"""meryl — k-mer counting and k-mer database operations.

Exposes:

- :func:`meryl_count` — count the k-mers of each reads file into its own
  database, ``<prefix>.<reads name>.meryl`` (the reads name without its
  ``.f*`` extension), as upstream.
- :func:`meryl_unionsum` — merge databases, summing the counts:
  ``<prefix>.unionsum.meryl``.
- :func:`meryl_histogram` — the k-mer count histogram of a database:
  ``<prefix>.hist``.

A meryl database is a directory of 129 files, so it is passed between tasks
as a ``<name>.meryl.tar.gz`` holding the ``<name>.meryl`` directory; the
tools unpack it under its name. ``threads`` follows the task's CPUs and
``memory`` (whole GiB) the task's memory: 90% of it for counting, as
upstream.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MERYL_IMAGE = "quay.io/biocontainers/meryl:1.4.1--h4ac6f70_0"

# Upstream's labels ask for 12 CPUs and 72 GB to count, 2 CPUs and 12 GB for
# the rest; these fit the test data and small genomes.
COUNT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")

# Unpack each database tarball (a <name>.meryl directory) into the working directory.
UNPACK = r"""
        unpack_dbs() {
            DBS=()
            for T in "$@"; do
                tar -xzf "$T"
                DBS+=("$(basename "$T" .tar.gz)")
            done
        }
"""

meryl_count_cmd = shell.create(
    name="meryl_count",
    image=MERYL_IMAGE,
    resources=COUNT_RESOURCES,
    inputs={"reads": list[File], "kvalue": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        MEMORY=$(( TASK_MEMORY_MB * 9 / 10 / 1024 ))
        for READ in "${{READS[@]}}"; do
            NAME=$(basename "$READ")
            DB="$P.${{NAME%.f*}}.meryl"
            meryl count k={{inputs.kvalue}} threads=$TASK_CPUS memory=$MEMORY $ARGS "$READ" output "$DB" || exit 1
            tar -czf "$OUT/$DB.tar.gz" "$DB"
        done
    """,
)

meryl_unionsum_cmd = shell.create(
    name="meryl_unionsum",
    image=MERYL_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"meryl_dbs": list[File], "kvalue": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.meryl_dbs}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {UNPACK}
        unpack_dbs "${{IN[@]}}"
        meryl union-sum k={{inputs.kvalue}} threads=$TASK_CPUS memory=$(( TASK_MEMORY_MB / 1024 )) $ARGS \
            output "$P.unionsum.meryl" "${{DBS[@]}}" || exit 1
        tar -czf {{outputs.results}}/"$P.unionsum.meryl.tar.gz" "$P.unionsum.meryl"
    """,
)

meryl_histogram_cmd = shell.create(
    name="meryl_histogram",
    image=MERYL_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"meryl_db": File, "kvalue": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.meryl_db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {UNPACK}
        unpack_dbs "${{IN[0]}}"
        meryl histogram k={{inputs.kvalue}} threads=$TASK_CPUS memory=$(( TASK_MEMORY_MB / 1024 )) $ARGS \
            "${{DBS[0]}}" > {{outputs.results}}/"$P.hist" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "meryl", meryl_count_cmd.as_task(), meryl_unionsum_cmd.as_task(), meryl_histogram_cmd.as_task()
)


async def meryl_count(reads: list[File], kvalue: int, prefix: str = "sample", args: str = "") -> list[File]:
    """One k-mer database per reads file: ``<prefix>.<reads name>.meryl.tar.gz``, in file-name order."""
    if not reads:
        raise ValueError("meryl count needs reads")
    results = await meryl_count_cmd(reads=reads, kvalue=kvalue, prefix=prefix, args=args)
    dbs = sorted([f async for f in results.walk() if f.name.endswith(".meryl.tar.gz")], key=lambda f: f.name)
    if len(dbs) != len(reads):
        raise FileNotFoundError(f"meryl count wrote {len(dbs)} databases for {len(reads)} reads files")
    return dbs


async def meryl_unionsum(meryl_dbs: list[File], kvalue: int, prefix: str = "sample", args: str = "") -> File:
    """The databases merged with summed counts: ``<prefix>.unionsum.meryl.tar.gz``."""
    results = await meryl_unionsum_cmd(meryl_dbs=meryl_dbs, kvalue=kvalue, prefix=prefix, args=args)
    db = await results.get_file(f"{prefix}.unionsum.meryl.tar.gz")
    if db is None:
        raise FileNotFoundError(f"meryl union-sum wrote no {prefix}.unionsum.meryl.tar.gz")
    return db


async def meryl_histogram(meryl_db: File, kvalue: int, prefix: str = "sample", args: str = "") -> File:
    """The k-mer count histogram of ``meryl_db``: ``<prefix>.hist``."""
    results = await meryl_histogram_cmd(meryl_db=meryl_db, kvalue=kvalue, prefix=prefix, args=args)
    hist = await results.get_file(f"{prefix}.hist")
    if hist is None:
        raise FileNotFoundError(f"meryl histogram wrote no {prefix}.hist")
    return hist
