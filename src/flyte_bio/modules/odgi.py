"""odgi — build, transform, lay out, draw and summarise variation graphs.

Exposes:

- :func:`odgi_build` — a GFA graph in odgi's format: ``<prefix>.og``.
- :func:`odgi_sort` — sort a graph's nodes (``args`` choose the
  algorithm, e.g. ``-p Ygs``): ``<prefix>.og``.
- :func:`odgi_unchop` — merge unitigs (nodes in a simple path) into single
  nodes: ``<prefix>.og``.
- :func:`odgi_squeeze` — merge several graphs into one: ``<prefix>.og``.
- :func:`odgi_view` — a graph as GFA: ``<prefix>.gfa``.
- :func:`odgi_layout` — a 2D layout of a graph: ``<prefix>.lay`` and, with
  ``tsv``, the same coordinates as ``<prefix>.tsv``.
- :func:`odgi_draw` — draw a graph's 2D layout: ``<prefix>.png``.
- :func:`odgi_viz` — a 1D visualisation of a graph's paths:
  ``<prefix>.png``.
- :func:`odgi_stats` — graph metrics: ``<prefix>.og.stats.tsv``, or
  ``<prefix>.og.stats.yaml`` when ``args`` has ``--yaml`` or
  ``--multiqc``.

Graphs may be given as ``.og`` or ``.gfa``; odgi reads either, telling them
apart by extension, so inputs are linked into the working directory under
their own names. As upstream, every tool gets the task's CPUs as threads.
Sorting and layout are multithreaded and stochastic, so their outputs
depend on the thread count. Options in ``args`` come last, except for
:func:`odgi_squeeze`, where they come first, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ODGI_IMAGE = "community.wave.seqera.io/library/odgi:0.9.0--48ccb29b570a8d65"

# Upstream labels build, layout, squeeze and unchop process_medium, sort
# process_high, and draw, view, viz and stats process_single.
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
HIGH_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Links the graph input $GRAPH_INPUT into a fresh working directory under its own name, as $GRAPH.
LINK_GRAPH = rf"""
        {TASK_LIMITS}
        set -e
        W=$(mktemp -d); cd "$W"
        GRAPH=$(basename "${{GRAPH_INPUT[0]}}")
        ln -s "${{GRAPH_INPUT[0]}}" "$GRAPH"
"""


def odgi_graph_cmd(name: str, resources: flyte.Resources, command: str):
    """A shell task running ``command`` on one graph, with inputs ``graph``, ``prefix`` and ``args``."""
    return shell.create(
        name=name,
        image=ODGI_IMAGE,
        resources=resources,
        inputs={"graph": File, "prefix": str, "args": str},
        defaults={"args": ""},
        outputs={"results": Dir},
        script=rf"""
        GRAPH_INPUT=({{inputs.graph}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {LINK_GRAPH}
        {command} || exit 1
        """,
    )


odgi_build_cmd = odgi_graph_cmd(
    "odgi_build", MEDIUM_RESOURCES, 'odgi build --threads $TASK_CPUS --gfa "$GRAPH" --out "$OUT/$P.og" $ARGS'
)
odgi_sort_cmd = odgi_graph_cmd(
    "odgi_sort", HIGH_RESOURCES, 'odgi sort --threads $TASK_CPUS --idx "$GRAPH" --out "$OUT/$P.og" $ARGS'
)
odgi_unchop_cmd = odgi_graph_cmd(
    "odgi_unchop", MEDIUM_RESOURCES, 'odgi unchop --threads $TASK_CPUS --idx "$GRAPH" --out "$OUT/$P.og" $ARGS'
)
odgi_view_cmd = odgi_graph_cmd(
    "odgi_view", SINGLE_RESOURCES, 'odgi view --threads $TASK_CPUS --idx "$GRAPH" --to-gfa $ARGS > "$OUT/$P.gfa"'
)
odgi_viz_cmd = odgi_graph_cmd(
    "odgi_viz", SINGLE_RESOURCES, 'odgi viz --threads $TASK_CPUS --idx "$GRAPH" --out "$OUT/$P.png" $ARGS'
)
# As upstream, --yaml or --multiqc in args make the output YAML.
odgi_stats_cmd = odgi_graph_cmd(
    "odgi_stats",
    SINGLE_RESOURCES,
    "SUFFIX=og.stats.tsv\n"
    '        case " $ARGS " in *--yaml*|*--multiqc*) SUFFIX=og.stats.yaml ;; esac\n'
    '        odgi stats --threads $TASK_CPUS --idx "$GRAPH" $ARGS > "$OUT/$P.$SUFFIX"',
)

odgi_layout_cmd = shell.create(
    name="odgi_layout",
    image=ODGI_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"graph": File, "prefix": str, "tsv": bool, "args": str},
    defaults={"tsv": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        GRAPH_INPUT=({{inputs.graph}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {LINK_GRAPH}
        TSV=()
        if [ {{inputs.tsv}} = true ]; then TSV=(--tsv "$OUT/$P.tsv"); fi
        odgi layout --threads $TASK_CPUS --idx "$GRAPH" --out "$OUT/$P.lay" "${{TSV[@]}}" $ARGS || exit 1
    """,
)

odgi_draw_cmd = shell.create(
    name="odgi_draw",
    image=ODGI_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"graph": File, "layout": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        GRAPH_INPUT=({{inputs.graph}})
        LAY=({{inputs.layout}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {LINK_GRAPH}
        ln -s "${{LAY[0]}}" "$(basename "${{LAY[0]}}")"
        odgi draw --threads $TASK_CPUS --idx "$GRAPH" --coords-in "$(basename "${{LAY[0]}}")" \
            --png "$OUT/$P.png" $ARGS || exit 1
    """,
)

# As upstream, the graphs are listed one name per line in a file for --input-graphs.
odgi_squeeze_cmd = shell.create(
    name="odgi_squeeze",
    image=ODGI_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"graphs": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        GRAPHS=({{inputs.graphs}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        : > files
        for f in "${{GRAPHS[@]}}"; do ln -s "$f" "$(basename "$f")"; basename "$f" >> files; done
        odgi squeeze $ARGS --threads $TASK_CPUS --input-graphs files -o "$OUT/$P.og" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "odgi",
    odgi_build_cmd.as_task(),
    odgi_sort_cmd.as_task(),
    odgi_unchop_cmd.as_task(),
    odgi_view_cmd.as_task(),
    odgi_viz_cmd.as_task(),
    odgi_stats_cmd.as_task(),
    odgi_layout_cmd.as_task(),
    odgi_draw_cmd.as_task(),
    odgi_squeeze_cmd.as_task(),
)


@dataclass
class OdgiLayoutResult:
    lay: File  # <prefix>.lay
    tsv: File | None  # <prefix>.tsv, with tsv=True
    results: Dir


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def odgi_build(graph: File, prefix: str = "sample", args: str = "") -> File:
    """``graph`` (GFA) in odgi's format: ``<prefix>.og``."""
    return await pick(await odgi_build_cmd(graph=graph, prefix=prefix, args=args), f"{prefix}.og", "odgi build")


async def odgi_sort(graph: File, prefix: str = "sample", args: str = "") -> File:
    """``graph`` with its nodes sorted: ``<prefix>.og``."""
    return await pick(await odgi_sort_cmd(graph=graph, prefix=prefix, args=args), f"{prefix}.og", "odgi sort")


async def odgi_unchop(graph: File, prefix: str = "sample", args: str = "") -> File:
    """``graph`` with its unitigs merged into single nodes: ``<prefix>.og``."""
    return await pick(await odgi_unchop_cmd(graph=graph, prefix=prefix, args=args), f"{prefix}.og", "odgi unchop")


async def odgi_squeeze(graphs: list[File], prefix: str = "sample", args: str = "") -> File:
    """``graphs`` merged into one graph: ``<prefix>.og``."""
    if not graphs:
        raise ValueError("odgi_squeeze: pass at least one graph")
    if any((g.name or g.path.rsplit("/", 1)[-1]) == f"{prefix}.og" for g in graphs):
        raise ValueError("odgi_squeeze: an input graph has the output's name; pass another prefix")
    return await pick(await odgi_squeeze_cmd(graphs=graphs, prefix=prefix, args=args), f"{prefix}.og", "odgi squeeze")


async def odgi_view(graph: File, prefix: str = "sample", args: str = "") -> File:
    """``graph`` as GFA: ``<prefix>.gfa``."""
    return await pick(await odgi_view_cmd(graph=graph, prefix=prefix, args=args), f"{prefix}.gfa", "odgi view")


async def odgi_viz(graph: File, prefix: str = "sample", args: str = "") -> File:
    """A 1D visualisation of ``graph``'s paths: ``<prefix>.png``."""
    return await pick(await odgi_viz_cmd(graph=graph, prefix=prefix, args=args), f"{prefix}.png", "odgi viz")


def odgi_stats_name(prefix: str, args: str) -> str:
    suffix = "og.stats.yaml" if "--yaml" in args or "--multiqc" in args else "og.stats.tsv"
    return f"{prefix}.{suffix}"


async def odgi_stats(graph: File, prefix: str = "sample", args: str = "") -> File:
    """Metrics of ``graph`` (``args`` choose them): ``<prefix>.og.stats.tsv`` or ``.og.stats.yaml``."""
    results = await odgi_stats_cmd(graph=graph, prefix=prefix, args=args)
    return await pick(results, odgi_stats_name(prefix, args), "odgi stats")


async def odgi_layout(graph: File, prefix: str = "sample", tsv: bool = False, args: str = "") -> OdgiLayoutResult:
    """A 2D layout of ``graph``: ``<prefix>.lay``, and ``<prefix>.tsv`` with ``tsv``."""
    results = await odgi_layout_cmd(graph=graph, prefix=prefix, tsv=tsv, args=args)
    return OdgiLayoutResult(
        lay=await pick(results, f"{prefix}.lay", "odgi layout"),
        tsv=await pick(results, f"{prefix}.tsv", "odgi layout") if tsv else None,
        results=results,
    )


async def odgi_draw(graph: File, layout: File, prefix: str = "sample", args: str = "") -> File:
    """Draw ``graph`` with its 2D ``layout`` (from :func:`odgi_layout`): ``<prefix>.png``."""
    results = await odgi_draw_cmd(graph=graph, layout=layout, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.png", "odgi draw")
