"""pbccs — PacBio's CCS: HiFi consensus reads from PacBio subreads.

Exposes :func:`pbccs`, which runs ``ccs`` on one chunk of a subreads BAM
(``--chunk <chunk_num>/<chunk_on>``, which needs the BAM's ``.pbi`` index)
and writes, as upstream:

- ``<prefix>.chunk<N>.bam`` and its index ``<prefix>.chunk<N>.bam.pbi``: the
  consensus reads;
- ``<prefix>.chunk<N>.report.txt`` and ``<prefix>.chunk<N>.report.json``:
  the ZMW yield report;
- ``<prefix>.chunk<N>.metrics.json.gz``: per-ZMW metrics.

As upstream, ccs gets the task's CPUs as ``-j``, and the BAM and its index
are linked into a working directory under their own names, so the BAM
header's ``@PG`` command line names them as upstream's does. ``args`` go to
ccs (the upstream pipeline uses ``--min-rq 0.9``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PBCCS_IMAGE = "quay.io/biocontainers/pbccs:6.4.0--h9ee0642_0"
# Upstream labels it process_low, but its test snapshot is of a 4-thread run: the
# output BAM's @PG command line records -j (reproduced with 4, not 1 or 2).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

pbccs_cmd = shell.create(
    name="pbccs",
    image=PBCCS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "pbi": File, "chunk_num": int, "chunk_on": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); PBI=({{inputs.pbi}})
        N={{inputs.chunk_num}}
        ON={{inputs.chunk_on}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$NAME"
        ln -s "${{PBI[0]}}" "$NAME.pbi"
        C="$P.chunk$N"
        ccs "$NAME" "$C.bam" \
            --report-file "$C.report.txt" \
            --report-json "$C.report.json" \
            --metrics-json "$C.metrics.json.gz" \
            --chunk "$N/$ON" \
            -j $TASK_CPUS \
            $ARGS || exit 1
        mv "$C.bam" "$C.bam.pbi" "$C.report.txt" "$C.report.json" "$C.metrics.json.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("pbccs", pbccs_cmd.as_task())


@dataclass
class PbccsResult:
    bam: File  # <prefix>.chunk<N>.bam
    pbi: File  # <prefix>.chunk<N>.bam.pbi
    report_txt: File  # <prefix>.chunk<N>.report.txt
    report_json: File  # <prefix>.chunk<N>.report.json
    metrics: File  # <prefix>.chunk<N>.metrics.json.gz
    results: Dir  # all of them


async def pbccs(
    bam: File, pbi: File, chunk_num: int = 1, chunk_on: int = 1, prefix: str = "sample", args: str = ""
) -> PbccsResult:
    """HiFi reads of chunk ``chunk_num`` of ``chunk_on`` of subreads ``bam`` (with its index ``pbi``)."""
    if not 1 <= chunk_num <= chunk_on:
        raise ValueError(f"chunk {chunk_num}/{chunk_on}: need 1 <= chunk_num <= chunk_on")
    results = await pbccs_cmd(bam=bam, pbi=pbi, chunk_num=chunk_num, chunk_on=chunk_on, prefix=prefix, args=args)
    stem = f"{prefix}.chunk{chunk_num}"
    files = {}
    for key, name in (
        ("bam", f"{stem}.bam"),
        ("pbi", f"{stem}.bam.pbi"),
        ("report_txt", f"{stem}.report.txt"),
        ("report_json", f"{stem}.report.json"),
        ("metrics", f"{stem}.metrics.json.gz"),
    ):
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"ccs wrote no {name}")
        files[key] = found
    return PbccsResult(results=results, **files)
