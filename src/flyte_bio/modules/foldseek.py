"""foldseek — protein structure search with Foldseek.

A Foldseek database is a set of files sharing a base name (``<name>``, its
``.dbtype``, ``.index``, ``.lookup`` and ``.source``, and the ``_ca``, ``_h``
and ``_ss`` databases). The functions here take and return a database as a Dir
holding them: a database built from a few structures holds a handful of files.
A large prebuilt database (PDB, AlphaFold DB) is tens of files of many GiB;
pass it as a Dir too only while it stays small enough to stage.

Exposes:

- :func:`foldseek_createdb` — a database from structure files (PDB or mmCIF,
  may be gzipped): a Dir of ``<prefix>``, ``<prefix>.{dbtype,index,lookup,source}``
  and ``<prefix>_{ca,h,ss}{,.dbtype,.index}``.
- :func:`foldseek_easysearch` — search query structures against a database:
  ``<prefix>.m8``, the hits in BLAST tabular format. As upstream, the database
  in the Dir is ``<db_name>`` if ``<db_name>.dbtype`` exists, else the single
  database found by its ``.lookup`` file.

Both run with ``--threads`` set to the task's CPUs. Temporary files go to a
scratch directory outside the outputs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FOLDSEEK_IMAGE = "quay.io/biocontainers/foldseek:10.941cd33--h5021889_1"
# Upstream: createdb process_single, easysearch process_medium. Searching a
# large database needs more memory and disk for its temporary files.
CREATEDB_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
EASYSEARCH_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")

foldseek_createdb_cmd = shell.create(
    name="foldseek_createdb",
    image=FOLDSEEK_IMAGE,
    resources=CREATEDB_RESOURCES,
    inputs={"structures": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"db": Dir},
    script=rf"""
        {TASK_LIMITS}
        STRUCTURES=({{inputs.structures}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.db}}
        W=$(mktemp -d); cd "$W"

        foldseek createdb "${{STRUCTURES[@]}}" "$OUT/$P" $ARGS --threads $TASK_CPUS || exit 1
    """,
)

foldseek_easysearch_cmd = shell.create(
    name="foldseek_easysearch",
    image=FOLDSEEK_IMAGE,
    resources=EASYSEARCH_RESOURCES,
    inputs={"query": list[File], "db": Dir, "db_name": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        QUERY=({{inputs.query}})
        DB_DIR={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}

        DB_NAME="$DB_DIR/"{{inputs.db_name}}
        if [ ! -f "$DB_NAME.dbtype" ]; then
            DB_NAME=$(find -L "$DB_DIR/" -maxdepth 1 -name '*.lookup' | sed 's/\.lookup$//')
            if [ "$(printf '%s\n' "$DB_NAME" | grep -c .)" -ne 1 ]; then
                echo "ERROR: $DB_DIR/{{inputs.db_name}} not found and expected exactly one .lookup file in $DB_DIR" >&2
                exit 1
            fi
        fi
        W=$(mktemp -d); cd "$W"

        foldseek easy-search "${{QUERY[@]}}" "$DB_NAME" "$OUT/$P.m8" tmp $ARGS --threads $TASK_CPUS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("foldseek", foldseek_createdb_cmd.as_task(), foldseek_easysearch_cmd.as_task())


async def foldseek_createdb(structures: File | list[File], prefix: str = "db", args: str = "") -> Dir:
    """A Foldseek database of ``structures``: a Dir of ``<prefix>*``."""
    if isinstance(structures, File):
        structures = [structures]
    return await foldseek_createdb_cmd(structures=structures, prefix=prefix, args=args)


async def foldseek_easysearch(
    query: File | list[File], db: Dir, db_name: str = "", prefix: str = "sample", args: str = ""
) -> File:
    """Search ``query`` structures against the database in ``db``: ``<prefix>.m8``.

    ``db_name`` is the database's base name in ``db`` (upstream's ``meta2.id``);
    if it isn't there, the single database in ``db`` is used.
    """
    if isinstance(query, File):
        query = [query]
    results = await foldseek_easysearch_cmd(query=query, db=db, db_name=db_name or "db", prefix=prefix, args=args)
    aln = await results.get_file(f"{prefix}.m8")
    if aln is None:
        raise FileNotFoundError(f"foldseek easy-search wrote no {prefix}.m8")
    return aln
