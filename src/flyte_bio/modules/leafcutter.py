"""leafcutter — intron clustering of splice junctions for differential splicing.

Exposes :func:`leafcutter_clusterregtools`, which runs ``leafcutter-cluster``
on junction files from ``regtools junctions extract`` (:mod:`regtools`) and
writes, as upstream, ``<prefix>_pooled`` (the pooled junctions),
``<prefix>_refined`` (the refined clusters), ``<prefix>_sortedlibs`` (the
per-sample sorted junction files), ``<prefix>_perind.counts.gz`` (per-sample
``count/cluster total`` ratios) and ``<prefix>_perind_numers.counts.gz``
(per-sample counts).

As upstream, the junction files are listed by name, sorted, and clustered in
the working directory (``-r ./``), where leafcutter also writes each sample's
``<junc file>.<prefix>.sorted.gz``. ``args`` takes leafcutter-cluster's
options (``-m``, ``-l``, ``-p``, ...).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

LEAFCUTTER_IMAGE = "quay.io/biocontainers/leafcutter:2.0.3--pyhd8ed1ab_0"
# Upstream: process_medium; leafcutter-cluster is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

leafcutter_clusterregtools_cmd = shell.create(
    name="leafcutter_clusterregtools",
    image=LEAFCUTTER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"junctions": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        JUNCTIONS=({inputs.junctions})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${JUNCTIONS[@]}"; do ln -s "$F" .; NAMES+=("$(basename "$F")"); done
        printf '%s\n' "${NAMES[@]}" | sort > juncfiles_list.txt
        leafcutter-cluster -j juncfiles_list.txt -r ./ -o "$P" $ARGS || exit 1
        mv "${P}_pooled" "${P}_refined" "${P}_sortedlibs" "${P}"_perind*.counts.gz "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("leafcutter", leafcutter_clusterregtools_cmd.as_task())


@dataclass
class LeafcutterClusters:
    pooled: File  # <prefix>_pooled
    refined: File  # <prefix>_refined
    sortedlibs: File  # <prefix>_sortedlibs
    counts: File  # <prefix>_perind.counts.gz
    numers: File  # <prefix>_perind_numers.counts.gz
    results: Dir


async def leafcutter_clusterregtools(junctions: list[File], prefix: str, args: str = "") -> LeafcutterClusters:
    """Cluster the introns of regtools ``junctions`` files (one per sample)."""
    results = await leafcutter_clusterregtools_cmd(junctions=junctions, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"leafcutter-cluster wrote no {name}")
        return found

    return LeafcutterClusters(
        pooled=await pick(f"{prefix}_pooled"),
        refined=await pick(f"{prefix}_refined"),
        sortedlibs=await pick(f"{prefix}_sortedlibs"),
        counts=await pick(f"{prefix}_perind.counts.gz"),
        numers=await pick(f"{prefix}_perind_numers.counts.gz"),
        results=results,
    )
