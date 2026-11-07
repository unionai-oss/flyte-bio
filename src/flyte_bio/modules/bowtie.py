"""bowtie — Bowtie (version 1) short-read alignment.

Exposes:

- :func:`bowtie1_build` — build a Bowtie index of a FASTA. The index files
  (``<prefix>.{1,2,3,4}.ebwt`` and ``<prefix>.rev.{1,2}.ebwt``) land in one
  output Dir.
- :func:`bowtie1_align` — align one sample's reads (one or two mates) as SAM,
  piped into ``samtools view`` for an unsorted BAM.

Outputs of :func:`bowtie1_align` keep upstream's names:

- the alignment: ``<prefix>.bam``;
- Bowtie's alignment summary (its stderr): ``<prefix>.out``;
- with ``save_unaligned``: ``<prefix>.unmapped.fastq.gz`` (single-end), or
  ``<prefix>.unmapped_1.fastq.gz`` and ``<prefix>.unmapped_2.fastq.gz``.

As upstream, Bowtie and samtools both use the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BOWTIE_IMAGE = "community.wave.seqera.io/library/bowtie_htslib_samtools:e1e242368ffcb5d3"
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")

bowtie1_build_cmd = shell.create(
    name="bowtie1_build",
    image=BOWTIE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        bowtie-build --threads $TASK_CPUS "${{FA[0]}}" {{outputs.index}}/{{inputs.prefix}} || exit 1
    """,
)

# As upstream: the index is found by its .3.ebwt file. Bowtie names the
# unaligned mates <prefix>.unmapped_{1,2}.fastq; they're compressed after.
# Bowtie's wrapper is Python and its errors may exit 255; any failure exits 1.
bowtie1_align_cmd = shell.create(
    name="bowtie1_align",
    image=BOWTIE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "prefix": str,
        "save_unaligned": bool,
        "args": str,
        "args2": str,
    },
    defaults={"save_unaligned": False, "args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}

        INDEX=$(find -L {{inputs.index}} -name "*.3.ebwt" | sed 's/\.3\.ebwt$//')
        if [ -z "$INDEX" ]; then echo "Bowtie index files not found" >&2; exit 1; fi

        UNALIGNED=(--un "$OUT/$P.unmapped.fastq")
        if [ {{inputs.save_unaligned}} != true ]; then UNALIGNED=(); fi
        if [ ${{#R2[@]}} -eq 0 ]; then READS=("${{R1[0]}}"); else READS=(-1 "${{R1[0]}}" -2 "${{R2[0]}}"); fi

        bowtie --threads $TASK_CPUS --sam -x "$INDEX" -q "${{UNALIGNED[@]}}" $ARGS "${{READS[@]}}" \
            2> "$OUT/$P.out" \
            | samtools view $ARGS2 -@ $TASK_CPUS -bS -o "$OUT/$P.bam" -
        STATUS=$?
        cat "$OUT/$P.out" >&2
        if [ $STATUS -ne 0 ]; then exit 1; fi

        for FQ in "$OUT/$P.unmapped.fastq" "$OUT/${{P}}.unmapped_1.fastq" "$OUT/${{P}}.unmapped_2.fastq"; do
            if [ -f "$FQ" ]; then gzip "$FQ"; fi
        done
        true
    """,
)

env = flyte.TaskEnvironment.from_task(
    "bowtie",
    bowtie1_build_cmd.as_task(),
    bowtie1_align_cmd.as_task(),
)


async def bowtie1_build(fasta: File, prefix: str) -> Dir:
    """Bowtie index of ``fasta``: a Dir of ``<prefix>.{1,2,3,4,rev.1,rev.2}.ebwt``."""
    return await bowtie1_build_cmd(fasta=fasta, prefix=prefix)


@dataclass
class Bowtie1AlignResult:
    bam: File  # <prefix>.bam
    log: File  # <prefix>.out
    unaligned: list[File]  # with save_unaligned: one file (single-end) or mates 1 and 2
    results: Dir  # everything the task wrote


async def bowtie1_align(
    reads_1: File,
    reads_2: File | None,
    index: Dir,
    prefix: str,
    save_unaligned: bool = False,
    args: str = "",
    args2: str = "",
) -> Bowtie1AlignResult:
    """Align one sample's reads with Bowtie; ``args`` go to Bowtie, ``args2`` to samtools.

    ``reads_2`` is None for single-end reads.
    """
    results = await bowtie1_align_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        prefix=prefix,
        save_unaligned=save_unaligned,
        args=args,
        args2=args2,
    )

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"bowtie align wrote no {name}")
        return file

    unaligned = []
    if save_unaligned:
        mates = ["_1", "_2"] if reads_2 is not None else [""]
        unaligned = [await pick(f"{prefix}.unmapped{mate}.fastq.gz") for mate in mates]
    return Bowtie1AlignResult(
        bam=await pick(f"{prefix}.bam"),
        log=await pick(f"{prefix}.out"),
        unaligned=unaligned,
        results=results,
    )
