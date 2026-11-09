"""isoseq — PacBio's Iso-Seq tools: full-length transcript reads from CCS reads.

Exposes :func:`isoseq_refine`, which runs ``isoseq refine`` on one sample's
primer-trimmed full-length reads (lima's output BAM): it trims poly(A) tails
(with ``--require-polya`` in ``args``) and removes concatemers, writing the
full-length non-chimeric (FLNC) reads and their reports as upstream:

- ``<prefix>.bam`` and its index ``<prefix>.bam.pbi``;
- ``<prefix>.consensusreadset.xml``: the PacBio dataset;
- ``<prefix>.filter_summary.report.json`` and ``<prefix>.report.csv``.

As upstream, isoseq gets the task's CPUs as ``-j``, and the inputs are linked
into a working directory under their own names, so the BAM header's ``@PG``
command line names them as upstream's does. ``args`` go to ``isoseq refine``.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ISOSEQ_IMAGE = "quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0"
# Upstream labels it process_low, but its test snapshot is of a 4-thread run: the
# output BAM's @PG command line records -j (reproduced with 4, not 1 or 2).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

isoseq_refine_cmd = shell.create(
    name="isoseq_refine",
    image=ISOSEQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "primers": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); PRIMERS=({{inputs.primers}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$NAME"
        PRIMERS_NAME=$(basename "${{PRIMERS[0]}}"); ln -s "${{PRIMERS[0]}}" "$PRIMERS_NAME"
        isoseq refine -j $TASK_CPUS $ARGS "$NAME" "$PRIMERS_NAME" "$P.bam" || exit 1
        mv "$P.bam" "$P.bam.pbi" "$P.consensusreadset.xml" "$P.filter_summary.report.json" "$P.report.csv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("isoseq", isoseq_refine_cmd.as_task())


@dataclass
class IsoseqRefineResult:
    bam: File  # <prefix>.bam: the FLNC reads
    pbi: File  # <prefix>.bam.pbi
    consensusreadset: File  # <prefix>.consensusreadset.xml
    summary: File  # <prefix>.filter_summary.report.json
    report: File  # <prefix>.report.csv
    results: Dir  # all of them


async def isoseq_refine(bam: File, primers: File, prefix: str = "sample", args: str = "") -> IsoseqRefineResult:
    """FLNC reads of ``bam`` (full-length reads, trimmed of ``primers`` by lima): ``<prefix>.bam``."""
    if (bam.name or os.path.basename(bam.path)) == f"{prefix}.bam":
        raise ValueError(f"input and output are both {prefix}.bam: set another prefix")
    results = await isoseq_refine_cmd(bam=bam, primers=primers, prefix=prefix, args=args)
    files = {}
    for key, name in (
        ("bam", f"{prefix}.bam"),
        ("pbi", f"{prefix}.bam.pbi"),
        ("consensusreadset", f"{prefix}.consensusreadset.xml"),
        ("summary", f"{prefix}.filter_summary.report.json"),
        ("report", f"{prefix}.report.csv"),
    ):
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"isoseq refine wrote no {name}")
        files[key] = found
    return IsoseqRefineResult(results=results, **files)
