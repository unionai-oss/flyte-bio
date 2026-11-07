"""pbjasmine — PacBio's Jasmine: 5mC CpG methylation calls from HiFi kinetics.

Exposes :func:`pbjasmine`, which runs ``jasmine`` on one PacBio BAM with
kinetics tags and writes the reads with ``MM``/``ML`` methylation tags to
``<prefix>.bam`` (upstream's prefix defaults to ``<id>_jasmine``).

As upstream, jasmine gets the task's CPUs as ``--num-threads`` and the input
is linked into a working directory under its own name, so the BAM header's
``@PG`` command line names it as upstream's does. ``args`` go to jasmine.
"""

import os

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PBJASMINE_IMAGE = "quay.io/biocontainers/pbjasmine:26.1.3--hd63eeec_0"
# Upstream labels it process_medium (6 CPUs). The output BAM's @PG command line
# records --num-threads.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

pbjasmine_cmd = shell.create(
    name="pbjasmine",
    image=PBJASMINE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$NAME"
        jasmine $ARGS --num-threads $TASK_CPUS "$NAME" "$P.bam" || exit 1
        mv "$P.bam" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("pbjasmine", pbjasmine_cmd.as_task())


async def pbjasmine(bam: File, prefix: str = "sample_jasmine", args: str = "") -> File:
    """``bam`` with 5mC CpG calls from its kinetics: ``<prefix>.bam``."""
    if (bam.name or os.path.basename(bam.path)) == f"{prefix}.bam":
        raise ValueError(f"input and output are both {prefix}.bam: set another prefix")
    results = await pbjasmine_cmd(bam=bam, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.bam")
    if out is None:
        raise FileNotFoundError(f"jasmine wrote no {prefix}.bam")
    return out
