"""pydamage — per-contig ancient-DNA damage estimates (pydamage analyze).

Exposes :func:`pydamage_analyze`, which fits pydamage's damage models to the
reads of an indexed BAM (the ``.bai`` is linked next to it) and returns
``<prefix>_pydamage_results.csv``: one row per reference (contig) with its
damage model parameters, the likelihood-ratio test's p-value and q-value, and
the prediction accuracy, as upstream (which renames pydamage's
``pydamage_results/pydamage_results.csv``). ``-p`` is the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PYDAMAGE_IMAGE = "community.wave.seqera.io/library/pydamage:1.0--1c195c7c48e87ebb"
# Upstream's process_medium is 6 CPUs and 36 GB; the test data needs far less.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

pydamage_analyze_cmd = shell.create(
    name="pydamage_analyze",
    image=PYDAMAGE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        BAI=({{inputs.bai}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" NUMBA_CACHE_DIR="$W/tmp" MPLCONFIGDIR="$W/tmp"
        NAME=$(basename "${{BAM[0]}}")
        ln -s "${{BAM[0]}}" "$NAME"
        ln -s "${{BAI[0]}}" "$NAME.bai"
        pydamage analyze $ARGS -p $TASK_CPUS "$NAME" || exit 1
        mv pydamage_results/pydamage_results.csv {{outputs.results}}/"${{P}}_pydamage_results.csv"
    """,
)

env = flyte.TaskEnvironment.from_task("pydamage", pydamage_analyze_cmd.as_task())


async def pydamage_analyze(bam: File, bai: File, prefix: str = "sample", args: str = "") -> File:
    """pydamage's per-reference damage estimates for ``bam``: ``<prefix>_pydamage_results.csv``."""
    results = await pydamage_analyze_cmd(bam=bam, bai=bai, prefix=prefix, args=args)
    name = f"{prefix}_pydamage_results.csv"
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"pydamage wrote no {name}")
    return out
