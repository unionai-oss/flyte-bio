"""diamond — protein sequence clustering with DIAMOND.

Exposes :func:`diamond_cluster`, which clusters the proteins of ``db`` (a
protein FASTA or a DIAMOND ``.dmnd`` database) with ``diamond cluster`` and
returns its cluster table, ``<prefix>.tsv``: one line per protein, its
cluster's representative then the protein.

As upstream, DIAMOND's memory limit (``-M``) is the task's memory and its
threads (``-p``) the task's CPUs. ``args`` passes the clustering options
(``--approx-id``, ``--member-cover``, sensitivity modes and so on).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DIAMOND_IMAGE = "community.wave.seqera.io/library/diamond:2.2.2--910ae987965d0d4d"
# Upstream's process_medium asks for 6 CPUs and 36 GB; a proteome of tens of
# thousands of proteins clusters in a few GiB. Raise with
# diamond_cluster_cmd.as_task().override(resources=...) for large databases.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

diamond_cluster_cmd = shell.create(
    name="diamond_cluster",
    image=DIAMOND_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        DB=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        diamond \
            cluster \
            $ARGS \
            -M "${{TASK_MEMORY_MB}}M" \
            -p $TASK_CPUS \
            -d "${{DB[0]}}" \
            -o "$OUT/$P.tsv"
    """,
)

env = flyte.TaskEnvironment.from_task("diamond", diamond_cluster_cmd.as_task())


async def diamond_cluster(db: File, prefix: str = "sample", args: str = "") -> File:
    """``diamond cluster`` of the proteins in ``db`` (FASTA or ``.dmnd``): ``<prefix>.tsv``."""
    results = await diamond_cluster_cmd(db=db, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.tsv")
    if found is None:
        raise FileNotFoundError(f"diamond cluster wrote no {prefix}.tsv")
    return found
