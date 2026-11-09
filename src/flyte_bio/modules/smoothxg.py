"""smoothxg — smooth a variation graph by realigning its local blocks with partial order alignment.

Exposes :func:`smoothxg`, which writes the smoothed graph
``<prefix>.smoothxg.gfa`` and, when ``args`` asks for one
(``--write-msa-in-maf-format <name>.maf``), the blocks' multiple alignment
as a MAF.

As upstream, smoothxg gets the task's CPUs as threads; options in ``args``
come last.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SMOOTHXG_IMAGE = "quay.io/biocontainers/smoothxg:0.8.0--h40c17d1_0"

# Upstream labels smoothxg process_high.
SMOOTHXG_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

smoothxg_cmd = shell.create(
    name="smoothxg",
    image=SMOOTHXG_IMAGE,
    resources=SMOOTHXG_RESOURCES,
    inputs={"gfa": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        GFA=({{inputs.gfa}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{GFA[0]}}")
        ln -s "${{GFA[0]}}" "$NAME"
        smoothxg --threads=$TASK_CPUS --gfa-in="$NAME" --smoothed-out="$P.smoothxg.gfa" $ARGS || exit 1
        mv "$P.smoothxg.gfa" "$OUT/"
        shopt -s nullglob
        for f in ./*.maf; do mv "$f" "$OUT/"; done
    """,
)

env = flyte.TaskEnvironment.from_task("smoothxg", smoothxg_cmd.as_task())


@dataclass
class SmoothxgResult:
    gfa: File  # <prefix>.smoothxg.gfa
    maf: File | None  # the MAF named in args' --write-msa-in-maf-format
    results: Dir


async def smoothxg(gfa: File, prefix: str = "sample", args: str = "") -> SmoothxgResult:
    """Smooth the graph ``gfa``: ``<prefix>.smoothxg.gfa``, plus a MAF when ``args`` asks for one."""
    results = await smoothxg_cmd(gfa=gfa, prefix=prefix, args=args)
    graph = await results.get_file(f"{prefix}.smoothxg.gfa")
    if graph is None:
        raise FileNotFoundError(f"smoothxg wrote no {prefix}.smoothxg.gfa")
    mafs = [f for f in await results.list_files() if f.path.endswith(".maf")]
    return SmoothxgResult(gfa=graph, maf=mafs[0] if mafs else None, results=results)
