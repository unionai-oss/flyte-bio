"""trimmomatic — adapter and quality trimming of Illumina reads.

Exposes :func:`trimmomatic`, which trims one sample's reads (single- or
paired-end) with the Trimmomatic steps in ``args`` (e.g.
``ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
MINLEN:36``; the bundled adapter files are found by name). Trimmomatic fails
when given no steps. ``quality_trim`` steps go before ``args``, as upstream's
second argument string does.

Outputs keep upstream's names:

- trimmed reads: ``<prefix>.SE.paired.trim.fastq.gz`` (single-end), or
  ``<prefix>.paired.trim_1.fastq.gz`` / ``<prefix>.paired.trim_2.fastq.gz``
  and the reads whose mate was dropped,
  ``<prefix>.unpaired.trim_1.fastq.gz`` / ``<prefix>.unpaired.trim_2.fastq.gz``;
- ``<prefix>_trim.log`` (the per-read trim log), ``<prefix>_out.log``
  (Trimmomatic's stderr) and ``<prefix>.summary``.

Trimmomatic's heap is 80% of the task's memory (upstream leaves the
launcher's default).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

TRIMMOMATIC_IMAGE = "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

trimmomatic_cmd = shell.create(
    name="trimmomatic",
    image=TRIMMOMATIC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "quality_trim": str, "args": str},
    defaults={"quality_trim": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        QUALITY_TRIM={{inputs.quality_trim}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        if [ ${{#R2[@]}} -gt 0 ]; then
            MODE=PE
            FILES=("${{R1[0]}}" "${{R2[0]}}"
                "$OUT/$P.paired.trim_1.fastq.gz" "$OUT/$P.unpaired.trim_1.fastq.gz"
                "$OUT/$P.paired.trim_2.fastq.gz" "$OUT/$P.unpaired.trim_2.fastq.gz")
        else
            MODE=SE
            FILES=("${{R1[0]}}" "$OUT/$P.SE.paired.trim.fastq.gz")
        fi
        STATUS=0
        trimmomatic -Xmx{JAVA_HEAP_MB}m \
            $MODE \
            -threads $TASK_CPUS \
            -trimlog "$OUT/${{P}}_trim.log" \
            -summary "$OUT/$P.summary" \
            "${{FILES[@]}}" \
            $QUALITY_TRIM \
            $ARGS 2> "$OUT/${{P}}_out.log" || STATUS=$?
        cat "$OUT/${{P}}_out.log" >&2
        if [ $STATUS -ne 0 ]; then exit 1; fi
    """,
)

env = flyte.TaskEnvironment.from_task("trimmomatic", trimmomatic_cmd.as_task())


@dataclass
class TrimmomaticResult:
    reads: list[File]  # the trimmed reads: one file, or mates 1 and 2
    unpaired_reads: list[File]  # paired-end: reads whose mate was dropped, from mates 1 and 2
    trim_log: File  # <prefix>_trim.log
    out_log: File  # <prefix>_out.log
    summary: File  # <prefix>.summary


async def trimmomatic(
    reads_1: File,
    reads_2: File | None = None,
    prefix: str = "sample",
    quality_trim: str = "",
    args: str = "",
) -> TrimmomaticResult:
    """Trim one sample's reads; ``reads_2`` is None for single-end."""
    results = await trimmomatic_cmd(
        reads_1=reads_1, reads_2=reads_2, prefix=prefix, quality_trim=quality_trim, args=args
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"trimmomatic wrote no {name}")
        return found

    if reads_2 is None:
        reads = [await pick(f"{prefix}.SE.paired.trim.fastq.gz")]
        unpaired = []
    else:
        reads = [await pick(f"{prefix}.paired.trim_{mate}.fastq.gz") for mate in (1, 2)]
        unpaired = [await pick(f"{prefix}.unpaired.trim_{mate}.fastq.gz") for mate in (1, 2)]
    return TrimmomaticResult(
        reads=reads,
        unpaired_reads=unpaired,
        trim_log=await pick(f"{prefix}_trim.log"),
        out_log=await pick(f"{prefix}_out.log"),
        summary=await pick(f"{prefix}.summary"),
    )
