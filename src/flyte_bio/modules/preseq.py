"""preseq — library complexity estimation (``preseq lc_extrap``).

Exposes :func:`preseq_lcextrap`, which extrapolates the number of distinct
reads a library would yield at greater sequencing depth, from one sample's
mapped reads (BED / mapped-read text by default; pass ``-B`` in ``args`` for
a BAM). Paired-end samples add ``-pe``, as upstream.

Outputs keep upstream's names: ``<prefix>.lc_extrap.txt`` (the
extrapolation curve) and ``<prefix>.command.log`` (preseq's stderr).

Upstream retries a failed run with ``-defects`` (skip the defect test), on
its second task attempt; here the retry happens within the same task, and
the log keeps both attempts' stderr.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PRESEQ_IMAGE = "quay.io/biocontainers/preseq:3.2.0--hdcf5f25_6"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

preseq_lcextrap_cmd = shell.create(
    name="preseq_lcextrap",
    image=PRESEQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "prefix": str, "single_end": bool, "args": str},
    defaults={"single_end": False, "args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        PE=()
        [ {inputs.single_end} = true ] || PE=(-pe)
        LOG="$OUT/$P.command.log"
        if ! preseq lc_extrap $ARGS "${PE[@]}" -output "$OUT/$P.lc_extrap.txt" "${BAM[0]}" 2> "$LOG"; then
            echo "preseq lc_extrap failed; retrying with -defects" >> "$LOG"
            preseq lc_extrap $ARGS -defects "${PE[@]}" -output "$OUT/$P.lc_extrap.txt" "${BAM[0]}" \
                2>> "$LOG" || { cat "$LOG" >&2; exit 1; }
        fi
        cat "$LOG" >&2
    """,
)

env = flyte.TaskEnvironment.from_task("preseq", preseq_lcextrap_cmd.as_task())


@dataclass
class PreseqLcextrapResult:
    lc_extrap: File  # <prefix>.lc_extrap.txt
    log: File  # <prefix>.command.log
    results: Dir


async def preseq_lcextrap(
    bam: File, prefix: str = "sample", single_end: bool = False, args: str = ""
) -> PreseqLcextrapResult:
    """``preseq lc_extrap`` of ``bam``'s mapped reads; ``-pe`` unless ``single_end``."""
    results = await preseq_lcextrap_cmd(bam=bam, prefix=prefix, single_end=single_end, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"preseq lc_extrap wrote no {name}")
        return found

    return PreseqLcextrapResult(
        lc_extrap=await pick(f"{prefix}.lc_extrap.txt"), log=await pick(f"{prefix}.command.log"), results=results
    )
