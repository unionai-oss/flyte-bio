"""gfaffix — collapse walk-preserving shared affixes (redundant prefixes and suffixes) of a variation graph.

Exposes :func:`gfaffix`, which writes the normalised graph
``<prefix>.gfaffix.gfa`` and the collapsed affixes ``<prefix>.affixes.txt``.
Options in ``args`` come first, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GFAFFIX_IMAGE = "quay.io/biocontainers/gfaffix:0.2.1--hc1c3326_0"

# Upstream labels gfaffix process_single.
GFAFFIX_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

gfaffix_cmd = shell.create(
    name="gfaffix",
    image=GFAFFIX_IMAGE,
    resources=GFAFFIX_RESOURCES,
    inputs={"gfa": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -e
        GFA=({inputs.gfa})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        # gfaffix exits 255 on errors, which the platform retries; report failure as 1.
        gfaffix $ARGS "${GFA[0]}" -o "$OUT/$P.gfaffix.gfa" -a "$OUT/$P.affixes.txt" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("gfaffix", gfaffix_cmd.as_task())


@dataclass
class GfaffixResult:
    gfa: File  # <prefix>.gfaffix.gfa
    affixes: File  # <prefix>.affixes.txt
    results: Dir


async def gfaffix(gfa: File, prefix: str = "sample", args: str = "") -> GfaffixResult:
    """Collapse the shared affixes of ``gfa``: ``<prefix>.gfaffix.gfa`` and ``<prefix>.affixes.txt``."""
    results = await gfaffix_cmd(gfa=gfa, prefix=prefix, args=args)
    graph = await results.get_file(f"{prefix}.gfaffix.gfa")
    affixes = await results.get_file(f"{prefix}.affixes.txt")
    if graph is None or affixes is None:
        raise FileNotFoundError(f"gfaffix wrote no {prefix}.gfaffix.gfa or {prefix}.affixes.txt")
    return GfaffixResult(gfa=graph, affixes=affixes, results=results)
