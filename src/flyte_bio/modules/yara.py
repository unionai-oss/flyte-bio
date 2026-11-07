"""yara — Yara exact read mapping.

Exposes:

- :func:`yara_index` — build a Yara index of a FASTA. The index files
  (``<fasta name>.{lf,rid,sa,txt}.*``, named after the whole FASTA file name,
  as upstream) land in one output Dir.
- :func:`yara_mapper` — map one sample's reads (one or two mates) and keep
  the mapped reads, sorted and indexed.

Outputs of :func:`yara_mapper` keep upstream's names:

- single-end: ``<prefix>.mapped.bam`` and ``<prefix>.mapped.bam.bai``;
- paired-end: each mate's mapped reads apart, ``<prefix>_1.mapped.bam`` and
  ``<prefix>_2.mapped.bam``, each with its ``.bai``.

Yara and samtools both use the task's CPUs. Yara writes reads in a
thread-dependent order, which ``samtools sort`` keeps for reads at the same
position, so the BAM's read order only reproduces at a fixed CPU count.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

YARA_IMAGE = "community.wave.seqera.io/library/yara:1.0.2--6d4e380b8f30498b"
YARA_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/samtools_yara:0fff9114c59346aa"
# Upstream's process_medium; enough for small genomes such as the HLA reference.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

yara_index_cmd = shell.create(
    name="yara_index",
    image=YARA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=r"""
        FA=({inputs.fasta})
        ARGS={inputs.args}
        cd "$(mktemp -d)"
        yara_indexer "${FA[0]}" -o {outputs.index}/"$(basename "${FA[0]}")" $ARGS || exit 1
    """,
)

# As upstream: the index prefix is the index file name without its last two
# extensions; single-end reads go straight through samtools, paired-end
# reads are split by mate.
yara_mapper_cmd = shell.create(
    name="yara_mapper",
    image=YARA_SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "index": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        INDEX=$(find -L {{inputs.index}} -name "*.lf.drp" | head -n 1 | sed 's/\.lf\.drp$//')
        if [ -z "$INDEX" ]; then echo "Yara index files not found" >&2; exit 1; fi

        if [ ${{#R2[@]}} -eq 0 ]; then
            yara_mapper $ARGS -t $TASK_CPUS -f bam "$INDEX" "${{R1[0]}}" \
                | samtools view -@ $TASK_CPUS -hb -F4 \
                | samtools sort -@ $TASK_CPUS > "$OUT/$P.mapped.bam" || exit 1
            samtools index -@ $TASK_CPUS "$OUT/$P.mapped.bam" || exit 1
        else
            yara_mapper $ARGS -t $TASK_CPUS -f bam "$INDEX" "${{R1[0]}}" "${{R2[0]}}" > output.bam || exit 1
            samtools view -@ $TASK_CPUS -hF 4 -f 0x40 -b output.bam \
                | samtools sort -@ $TASK_CPUS > "$OUT/${{P}}_1.mapped.bam" || exit 1
            samtools view -@ $TASK_CPUS -hF 4 -f 0x80 -b output.bam \
                | samtools sort -@ $TASK_CPUS > "$OUT/${{P}}_2.mapped.bam" || exit 1
            samtools index -@ $TASK_CPUS "$OUT/${{P}}_1.mapped.bam" || exit 1
            samtools index -@ $TASK_CPUS "$OUT/${{P}}_2.mapped.bam" || exit 1
        fi
    """,
)

yara_index_env = flyte.TaskEnvironment.from_task("yara_index", yara_index_cmd.as_task())
yara_mapper_env = flyte.TaskEnvironment.from_task("yara_mapper", yara_mapper_cmd.as_task())
env = flyte.TaskEnvironment(name="yara", depends_on=[yara_index_env, yara_mapper_env])


async def yara_index(fasta: File, args: str = "") -> Dir:
    """Yara index of ``fasta``: a Dir of ``<fasta file name>.{lf,rid,sa,txt}.*``."""
    return await yara_index_cmd(fasta=fasta, args=args)


@dataclass
class YaraMapperResult:
    bams: list[File]  # <prefix>.mapped.bam, or <prefix>_1.mapped.bam and <prefix>_2.mapped.bam
    bais: list[File]  # each BAM's .bai
    results: Dir  # everything the task wrote


async def yara_mapper(reads_1: File, reads_2: File | None, index: Dir, prefix: str, args: str = "") -> YaraMapperResult:
    """Map one sample's reads with Yara, keeping the mapped reads; ``reads_2`` is None for single-end."""
    results = await yara_mapper_cmd(reads_1=reads_1, reads_2=reads_2, index=index, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"yara mapper wrote no {name}")
        return file

    names = [f"{prefix}.mapped.bam"] if reads_2 is None else [f"{prefix}_1.mapped.bam", f"{prefix}_2.mapped.bam"]
    return YaraMapperResult(
        bams=[await pick(name) for name in names],
        bais=[await pick(f"{name}.bai") for name in names],
        results=results,
    )
