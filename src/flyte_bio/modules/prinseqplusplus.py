"""prinseqplusplus — quality filtering of reads with PRINSEQ++.

Exposes :func:`prinseqplusplus`, which filters one sample's reads
(single- or paired-end) by the criteria in ``args`` (e.g. ``-lc_entropy=0.8``
for low complexity) and returns, under upstream's names:

- the reads that pass: ``<prefix>_good_out.fastq.gz``, or
  ``<prefix>_good_out_R1.fastq.gz`` / ``<prefix>_good_out_R2.fastq.gz``;
- paired-end reads whose mate failed: ``<prefix>_single_out_R{1,2}.fastq.gz``;
- the reads that fail: ``<prefix>_bad_out.fastq.gz`` or
  ``<prefix>_bad_out_R{1,2}.fastq.gz``;
- PRINSEQ++'s summary, ``<prefix>.log`` (its verbose standard output).

As upstream, ``-threads`` is the task's CPUs and outputs are gzipped.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PRINSEQPLUSPLUS_IMAGE = "quay.io/biocontainers/prinseq-plus-plus:1.2.3--hc90279e_1"
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

prinseqplusplus_cmd = shell.create(
    name="prinseqplusplus",
    image=PRINSEQPLUSPLUS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        INPUT=(-fastq "${{R1[0]}}")
        if [ ${{#R2[@]}} -gt 0 ]; then INPUT+=(-fastq2 "${{R2[0]}}"); fi
        W=$(mktemp -d); cd "$W"
        prinseq++ \
            -threads $TASK_CPUS \
            "${{INPUT[@]}}" \
            -out_name "$P" \
            -out_gz \
            -VERBOSE 1 \
            $ARGS \
            | tee "$P.log" \
            || exit 1
        shopt -s nullglob
        mv ./*_good_out*.fastq.gz ./*_single_out*.fastq.gz ./*_bad_out*.fastq.gz "$P.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("prinseqplusplus", prinseqplusplus_cmd.as_task())


@dataclass
class PrinseqplusplusResult:
    good_reads: list[File]  # one file, or R1 and R2
    single_reads: list[File]  # paired-end reads whose mate failed (R1, R2), if any
    bad_reads: list[File]  # one file, or R1 and R2, if any
    log: File


async def prinseqplusplus(
    reads_1: File,
    reads_2: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> PrinseqplusplusResult:
    """Filter one sample's reads; ``reads_2`` is None for single-end reads."""
    results = await prinseqplusplus_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)
    files = sorted(await results.list_files(), key=lambda f: f.name)

    def matching(kind: str) -> list[File]:
        return [f for f in files if f.name.startswith(f"{prefix}_{kind}_out") and f.name.endswith(".fastq.gz")]

    logs = [f for f in files if f.name == f"{prefix}.log"]
    good = matching("good")
    if not good or not logs:
        raise FileNotFoundError(f"prinseq++ wrote no {prefix}_good_out*.fastq.gz or {prefix}.log")
    return PrinseqplusplusResult(
        good_reads=good, single_reads=matching("single"), bad_reads=matching("bad"), log=logs[0]
    )
