"""clusty — clustering of large sparse distance or similarity graphs.

Exposes :func:`clusty`, which clusters the pairs in a table of distances or
similarities (e.g. ANI between genomes) with Clusty and returns the cluster
assignments, ``<prefix>.tsv``. ``args`` names the table's columns and the
algorithm's settings (e.g. ``--id-cols id1 id2 --distance-col ani
--similarity --min ani 0.70``). With ``objects``, a file listing every object
to cluster, objects with no pairs become singleton clusters.

As upstream, ``-t`` is the task's CPUs, and ``prefix`` mustn't make the
output's name the same as the table's.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CLUSTY_IMAGE = "quay.io/biocontainers/clusty:1.2.2--h9ee0642_0"
# Upstream: process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

clusty_cmd = shell.create(
    name="clusty",
    image=CLUSTY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"distances": File, "objects": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        DIST=({{inputs.distances}})
        shopt -s nullglob; OBJECTS=({{inputs.objects}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        if [ "$(basename "${{DIST[0]}}")" = "$P.tsv" ]; then
            echo "Input and output names are the same, set prefix to disambiguate!" >&2; exit 1
        fi
        OBJECTS_ARG=()
        if [ ${{#OBJECTS[@]}} -gt 0 ]; then OBJECTS_ARG=(--objects-file "${{OBJECTS[0]}}"); fi
        W=$(mktemp -d); cd "$W"
        clusty \
            $ARGS \
            -t $TASK_CPUS \
            "${{OBJECTS_ARG[@]}}" \
            "${{DIST[0]}}" \
            "$P.tsv" \
            || exit 1
        mv "$P.tsv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("clusty", clusty_cmd.as_task())


async def clusty(distances: File, objects: File | None = None, prefix: str = "sample", args: str = "") -> File:
    """Cluster the pairs in ``distances`` (and ``objects``' singletons): the assignments ``<prefix>.tsv``."""
    results = await clusty_cmd(distances=distances, objects=objects, prefix=prefix, args=args)
    assignments = await results.get_file(f"{prefix}.tsv")
    if assignments is None:
        raise FileNotFoundError(f"clusty wrote no {prefix}.tsv")
    return assignments
