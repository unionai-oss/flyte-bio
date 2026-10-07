"""ucsc — UCSC genome browser utilities.

Exposes :func:`bedclip`, which drops bedGraph records outside the
chromosome bounds in a chrom.sizes file, and :func:`bedgraphtobigwig`, which
converts a (clipped, sorted) bedGraph to bigWig. Outputs keep upstream's
names: ``<prefix>.bedGraph`` and ``<prefix>.bigWig``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URIs. Upstream's main.nf names the bare
# `biocontainers/...` (Docker Hub) — use the quay.io mirrors.
BEDCLIP_IMAGE = "quay.io/biocontainers/ucsc-bedclip:377--h0b8a92a_2"
BEDGRAPHTOBIGWIG_IMAGE = "quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

bedclip_cmd = shell.create(
    name="ucsc_bedclip",
    image=BEDCLIP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bedgraph": File, "sizes": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        bedClip {inputs.bedgraph} {inputs.sizes} {outputs.results}/{inputs.prefix}.bedGraph
    """,
)

bedgraphtobigwig_cmd = shell.create(
    name="ucsc_bedgraphtobigwig",
    image=BEDGRAPHTOBIGWIG_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bedgraph": File, "sizes": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        bedGraphToBigWig {inputs.bedgraph} {inputs.sizes} {outputs.results}/{inputs.prefix}.bigWig
    """,
)


# One image per tool, so one env each; `env` aggregates them.
bedclip_env = flyte.TaskEnvironment.from_task("ucsc_bedclip", bedclip_cmd.as_task())
bedgraphtobigwig_env = flyte.TaskEnvironment.from_task("ucsc_bedgraphtobigwig", bedgraphtobigwig_cmd.as_task())
env = flyte.TaskEnvironment(name="ucsc", depends_on=[bedclip_env, bedgraphtobigwig_env])


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"no {name}")
    return found


async def bedclip(bedgraph: File, sizes: File, prefix: str) -> File:
    """Clip ``bedgraph`` to the chromosome sizes; returns ``<prefix>.bedGraph``."""
    return await pick(await bedclip_cmd(bedgraph=bedgraph, sizes=sizes, prefix=prefix), f"{prefix}.bedGraph")


async def bedgraphtobigwig(bedgraph: File, sizes: File, prefix: str) -> File:
    """Convert ``bedgraph`` to ``<prefix>.bigWig``."""
    return await pick(
        await bedgraphtobigwig_cmd(bedgraph=bedgraph, sizes=sizes, prefix=prefix), f"{prefix}.bigWig"
    )
