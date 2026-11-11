"""adapterremoval — adapter trimming, quality trimming and read merging (AdapterRemoval 2).

Exposes :func:`adapterremoval`, which trims one sample's reads (single- or
paired-end) and, with ``--collapse`` in ``args``, merges overlapping pairs.
Adapters can come from an adapter list (``--adapter-list``). As upstream the
random seed is fixed (``--seed 42``), the outputs are gzipped, and
AdapterRemoval's ``.gz`` names are renamed to ``.fastq.gz``:

- ``<prefix>.truncated.fastq.gz``: single-end reads, or paired reads whose
  mate was discarded;
- ``<prefix>.discarded.fastq.gz``: discarded reads;
- ``<prefix>.pair1.truncated.fastq.gz`` / ``<prefix>.pair2.truncated.fastq.gz``:
  paired-end reads;
- ``<prefix>.collapsed.fastq.gz`` / ``<prefix>.collapsed.truncated.fastq.gz``:
  merged pairs (with ``--collapse``), untrimmed / quality-trimmed;
- ``<prefix>.paired.fastq.gz``: interleaved pairs (with ``--interleaved-output``);
- ``<prefix>.settings``: the run's settings and statistics.

``--threads`` is the task's CPUs; other options go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ADAPTERREMOVAL_IMAGE = "community.wave.seqera.io/library/adapterremoval:2.3.4--66c0d06f6b0408b0"
# Upstream labels it process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

adapterremoval_cmd = shell.create(
    name="adapterremoval",
    image=ADAPTERREMOVAL_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "adapter_list": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); LIST=({{inputs.adapter_list}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        INPUT=(--file1 "${{R1[0]}}")
        if [ ${{#R2[@]}} -gt 0 ]; then INPUT+=(--file2 "${{R2[0]}}"); fi
        ADAPTERS=()
        if [ ${{#LIST[@]}} -gt 0 ]; then ADAPTERS=(--adapter-list "${{LIST[0]}}"); fi
        AdapterRemoval \
            "${{INPUT[@]}}" \
            $ARGS \
            "${{ADAPTERS[@]}}" \
            --basename "$P" \
            --threads $TASK_CPUS \
            --seed 42 \
            --gzip \
            || exit 1
        for KIND in truncated discarded pair1.truncated pair2.truncated collapsed collapsed.truncated paired; do
            if [ -f "$P.$KIND.gz" ]; then mv "$P.$KIND.gz" "$P.$KIND.fastq.gz"; fi
        done
        shopt -s nullglob
        mv ./*.fastq.gz ./*.settings "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("adapterremoval", adapterremoval_cmd.as_task())


@dataclass
class AdapterRemovalResult:
    settings: File  # <prefix>.settings
    singles_truncated: File | None  # <prefix>.truncated.fastq.gz
    discarded: File | None  # <prefix>.discarded.fastq.gz
    paired_truncated: list[File]  # <prefix>.pair{1,2}.truncated.fastq.gz, paired-end
    collapsed: File | None  # <prefix>.collapsed.fastq.gz, with --collapse
    collapsed_truncated: File | None  # <prefix>.collapsed.truncated.fastq.gz, with --collapse
    paired_interleaved: File | None  # <prefix>.paired.fastq.gz, with --interleaved-output
    results: Dir  # everything AdapterRemoval wrote


async def adapterremoval(
    reads_1: File,
    reads_2: File | None = None,
    adapter_list: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> AdapterRemovalResult:
    """Trim one sample's reads (``reads_2`` None for single-end) with AdapterRemoval."""
    results = await adapterremoval_cmd(
        reads_1=reads_1, reads_2=reads_2, adapter_list=adapter_list, prefix=prefix, args=args
    )
    settings = await results.get_file(f"{prefix}.settings")
    if settings is None:
        raise FileNotFoundError(f"AdapterRemoval wrote no {prefix}.settings")

    async def pick(kind: str) -> File | None:
        return await results.get_file(f"{prefix}.{kind}.fastq.gz")

    pairs = [await pick("pair1.truncated"), await pick("pair2.truncated")]
    return AdapterRemovalResult(
        settings=settings,
        singles_truncated=await pick("truncated"),
        discarded=await pick("discarded"),
        paired_truncated=[f for f in pairs if f is not None],
        collapsed=await pick("collapsed"),
        collapsed_truncated=await pick("collapsed.truncated"),
        paired_interleaved=await pick("paired"),
        results=results,
    )
