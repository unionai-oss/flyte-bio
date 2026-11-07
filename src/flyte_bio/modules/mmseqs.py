"""mmseqs — sequence search and clustering with MMseqs2.

An MMseqs2 database is a set of files sharing a base name (``<name>``, its
``.index``, ``.dbtype``, ``.lookup``, a ``_h`` header database, ...). The
functions here take and return databases as a Dir holding them. As upstream,
a database in an input Dir is found by its files: the longest common prefix of
the files matching a pattern (default ``*.dbtype``), extension removed. So a
Dir may hold several databases with a common prefix, e.g. ``db``, ``db_h`` and
``db_clu``, and the pattern picks one (``*_clu.dbtype``) when the prefix
alone wouldn't. Databases hold a few files each (a clustering result one per
thread), so they're passed as a Dir rather than an archive.

Exposes:

- :func:`mmseqs_createdb` — a database from a FASTA (nucleotide or protein,
  may be gzipped): a Dir of ``<prefix>``, ``<prefix>.{dbtype,index,lookup,source}``
  and ``<prefix>_h{,.dbtype,.index}``.
- :func:`mmseqs_cluster` — cluster a database (cascaded clustering): a Dir
  of the cluster database ``<prefix>.{0,1,...,dbtype,index}``.
- :func:`mmseqs_linclust` — the same with linear-time clustering.
- :func:`mmseqs_createtsv` — a result database (clustering, search,
  taxonomy, ...) as ``<prefix>.tsv``, against its query (and target) database.
- :func:`mmseqs_databases` — download and build a public database
  (``SILVA``, ``UniRef50``, ...; see ``mmseqs databases -h``): a Dir of
  ``database*``. Needs network access; large databases need the task's disk
  and memory raised.
- :func:`mmseqs_taxonomy` — assign taxonomy to a query database against a
  database with taxonomy (e.g. from :func:`mmseqs_databases`): a Dir of the
  taxonomy result database ``<prefix>.*``.
- :func:`mmseqs_createindex` — precompute the k-mer index of a database for
  faster searches: the input Dir with ``<name>.idx{,.dbtype,.index}`` added.
- :func:`mmseqs_tsv2exprofiledb` — an expandable profile database (as
  ColabFold's ``uniref30``/``colabfold_envdb`` downloads ship) from its
  ``<name>.tsv``, ``<name>_h.tsv``, ``<name>_seq.tsv`` and ``<name>_aln.tsv``:
  the input Dir with ``<name>_db``, ``<name>_db_h``, ``<name>_db_seq``,
  ``<name>_db_seq_h`` and ``<name>_db_aln`` databases added.

As upstream, createindex and tsv2exprofiledb return their input directory
with the new files beside the database, which is copied into the output (so
the task's disk needs room for the database twice). An index or profile
database is a few large files, not many, so it stays a Dir.

Every tool but createdb and tsv2exprofiledb runs with ``--threads`` set to the task's CPUs. A
clustering result has one data file per thread, so its file list follows the
CPU count. Temporary files go to a scratch directory outside the outputs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MMSEQS_IMAGE = "community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139"
# mmseqs databases also needs wget.
MMSEQS_DOWNLOAD_IMAGE = "community.wave.seqera.io/library/mmseqs2_wget:aa683a2c5355899d"
# Upstream: createdb process_low, createtsv process_single, the rest process_medium/high.
# Clustering or searching large databases needs more memory (MMseqs2 splits its
# prefilter to fit) and disk for its temporary files.
CREATEDB_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
CLUSTER_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")
CREATETSV_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
TAXONOMY_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")
DATABASES_RESOURCES = flyte.Resources(cpu=2, memory="8Gi", disk="40Gi")
# Upstream: process_high (+ process_high_memory for createindex). Real profile
# databases need far more memory and disk; raise them with an override.
CREATEINDEX_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="40Gi")
TSV2EXPROFILEDB_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="40Gi")
DEFAULT_PATTERN = "*.dbtype"

# The base path of the database in directory $1 whose files match $2, as
# upstream finds it: the matching files, extension removed, reduced to their
# longest common prefix.
DATABASE_BASE = r"""
        database_base() {
            local BASE
            BASE=$(find -L "$1/" -maxdepth 1 -name "$2" | sed 's/\.[^.]*$//' | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')
            if [ -z "$BASE" ]; then
                echo "no database files matching $2 in $1" >&2
                return 1
            fi
            echo "$BASE"
        }
""".strip()

mmseqs_createdb_cmd = shell.create(
    name="mmseqs_createdb",
    image=MMSEQS_IMAGE,
    resources=CREATEDB_RESOURCES,
    inputs={"sequence": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SEQ=({inputs.sequence})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        # The database's .source records the FASTA's name, so it's read under
        # its own name (decompressed, as upstream, if gzipped).
        NAME=$(basename "${SEQ[0]}")
        if [ "${NAME##*.}" = gz ]; then
            NAME=${NAME%.gz}
            gzip -c -d "${SEQ[0]}" > "$NAME"
        else
            ln -s "${SEQ[0]}" "$NAME"
        fi
        mmseqs createdb "$NAME" "$OUT/$P" $ARGS
    """,
)


def clustering_cmd(tool: str):
    """``mmseqs <tool>`` (cluster or linclust) of the database in ``db``."""
    return shell.create(
        name=f"mmseqs_{tool}",
        image=MMSEQS_IMAGE,
        resources=CLUSTER_RESOURCES,
        inputs={"db": Dir, "db_pattern": str, "prefix": str, "args": str},
        defaults={"db_pattern": DEFAULT_PATTERN, "args": ""},
        outputs={"results": Dir},
        script=rf"""
            {TASK_LIMITS}
            {DATABASE_BASE}
            P={{inputs.prefix}}
            ARGS={{inputs.args}}
            OUT={{outputs.results}}
            DB=$(database_base {{inputs.db}} {{inputs.db_pattern}}) || exit 1
            W=$(mktemp -d); cd "$W"

            mmseqs {tool} "$DB" "$OUT/$P" tmp1 $ARGS --threads $TASK_CPUS
        """,
    )


mmseqs_cluster_cmd = clustering_cmd("cluster")
mmseqs_linclust_cmd = clustering_cmd("linclust")

mmseqs_createtsv_cmd = shell.create(
    name="mmseqs_createtsv",
    image=MMSEQS_IMAGE,
    resources=CREATETSV_RESOURCES,
    inputs={
        "result": Dir,
        "query": Dir | None,
        "target": Dir,
        "result_pattern": str,
        "query_pattern": str,
        "target_pattern": str,
        "prefix": str,
        "args": str,
    },
    defaults={
        "result_pattern": DEFAULT_PATTERN,
        "query_pattern": DEFAULT_PATTERN,
        "target_pattern": DEFAULT_PATTERN,
        "args": "",
    },
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {DATABASE_BASE}
        QUERY_DIR={{inputs.query}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}

        # As upstream, without a query database the target takes its place.
        DBS=()
        if [ -d "$QUERY_DIR" ]; then
            QUERY=$(database_base "$QUERY_DIR" {{inputs.query_pattern}}) || exit 1
            DBS+=("$QUERY")
        fi
        TARGET=$(database_base {{inputs.target}} {{inputs.target_pattern}}) || exit 1
        RESULT=$(database_base {{inputs.result}} {{inputs.result_pattern}}) || exit 1
        DBS+=("$TARGET" "$RESULT")
        W=$(mktemp -d); cd "$W"

        mmseqs createtsv "${{DBS[@]}}" "$OUT/$P.tsv" $ARGS --threads $TASK_CPUS
    """,
)

mmseqs_taxonomy_cmd = shell.create(
    name="mmseqs_taxonomy",
    image=MMSEQS_IMAGE,
    resources=TAXONOMY_RESOURCES,
    inputs={
        "query": Dir,
        "target": Dir,
        "query_pattern": str,
        "target_pattern": str,
        "prefix": str,
        "args": str,
    },
    defaults={"query_pattern": DEFAULT_PATTERN, "target_pattern": DEFAULT_PATTERN, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {DATABASE_BASE}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        QUERY=$(database_base {{inputs.query}} {{inputs.query_pattern}}) || exit 1
        TARGET=$(database_base {{inputs.target}} {{inputs.target_pattern}}) || exit 1
        W=$(mktemp -d); cd "$W"

        mmseqs taxonomy "$QUERY" "$TARGET" "$OUT/$P" tmp1 $ARGS --threads $TASK_CPUS
    """,
)

mmseqs_databases_cmd = shell.create(
    name="mmseqs_databases",
    image=MMSEQS_DOWNLOAD_IMAGE,
    resources=DATABASES_RESOURCES,
    inputs={"database": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        NAME={{inputs.database}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        mmseqs databases "$NAME" "$OUT/database" tmp/ --threads $TASK_CPUS $ARGS
    """,
)

mmseqs_createindex_cmd = shell.create(
    name="mmseqs_createindex",
    image=MMSEQS_IMAGE,
    resources=CREATEINDEX_RESOURCES,
    inputs={"db": Dir, "db_pattern": str, "args": str},
    defaults={"db_pattern": DEFAULT_PATTERN, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {DATABASE_BASE}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        # The index goes next to the database, so the database is copied into the output first.
        cp -rL {{inputs.db}}/. "$OUT"/ || exit 1
        DB=$(database_base "$OUT" {{inputs.db_pattern}}) || exit 1
        W=$(mktemp -d); cd "$W"

        # As upstream, the split memory limit is 80% of the task's memory, in whole GiB.
        mmseqs createindex "$DB" tmp1 $ARGS --threads $TASK_CPUS \
            --split-memory-limit $(( TASK_MEMORY_MB / 1024 * 8 / 10 ))G || exit 1
    """,
)

mmseqs_tsv2exprofiledb_cmd = shell.create(
    name="mmseqs_tsv2exprofiledb",
    image=MMSEQS_IMAGE,
    resources=TSV2EXPROFILEDB_RESOURCES,
    inputs={"database": Dir, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ARGS={inputs.args}
        OUT={outputs.results}
        cp -rL {inputs.database}/. "$OUT"/ || exit 1
        # As upstream, the database is named by its <name>_seq.tsv.
        DB=$(find -L "$OUT/" -name "*_seq.tsv" | sed 's/_seq\.tsv$//')
        if [ -z "$DB" ]; then
            echo "no *_seq.tsv in the database directory" >&2
            exit 1
        fi
        W=$(mktemp -d); cd "$W"

        mmseqs tsv2exprofiledb "$DB" "${DB}_db" $ARGS || exit 1
    """,
)

# One image per environment: mmseqs databases needs wget too.
mmseqs_env = flyte.TaskEnvironment.from_task(
    "mmseqs_tools",
    mmseqs_createdb_cmd.as_task(),
    mmseqs_cluster_cmd.as_task(),
    mmseqs_linclust_cmd.as_task(),
    mmseqs_createtsv_cmd.as_task(),
    mmseqs_taxonomy_cmd.as_task(),
    mmseqs_createindex_cmd.as_task(),
    mmseqs_tsv2exprofiledb_cmd.as_task(),
)
mmseqs_download_env = flyte.TaskEnvironment.from_task("mmseqs_download", mmseqs_databases_cmd.as_task())
env = flyte.TaskEnvironment(name="mmseqs", depends_on=[mmseqs_env, mmseqs_download_env])


async def mmseqs_createdb(sequence: File, prefix: str = "sample", args: str = "") -> Dir:
    """A database of ``sequence`` (may be gzipped): ``<prefix>``, ``<prefix>.*`` and ``<prefix>_h*``."""
    return await mmseqs_createdb_cmd(sequence=sequence, prefix=prefix, args=args)


async def mmseqs_cluster(db: Dir, prefix: str = "sample", db_pattern: str = DEFAULT_PATTERN, args: str = "") -> Dir:
    """Cluster the database in ``db`` (found by ``db_pattern``): a Dir of ``<prefix>.*``."""
    return await mmseqs_cluster_cmd(db=db, db_pattern=db_pattern, prefix=prefix, args=args)


async def mmseqs_linclust(db: Dir, prefix: str = "sample", db_pattern: str = DEFAULT_PATTERN, args: str = "") -> Dir:
    """Cluster the database in ``db`` in linear time (found by ``db_pattern``): a Dir of ``<prefix>.*``."""
    return await mmseqs_linclust_cmd(db=db, db_pattern=db_pattern, prefix=prefix, args=args)


async def mmseqs_createtsv(
    result: Dir,
    query: Dir | None,
    target: Dir,
    prefix: str = "sample",
    result_pattern: str = DEFAULT_PATTERN,
    query_pattern: str = DEFAULT_PATTERN,
    target_pattern: str = DEFAULT_PATTERN,
    args: str = "",
) -> File:
    """The result database in ``result`` as ``<prefix>.tsv``; without ``query``, ``target`` is the query."""
    results = await mmseqs_createtsv_cmd(
        result=result,
        query=query,
        target=target,
        result_pattern=result_pattern,
        query_pattern=query_pattern,
        target_pattern=target_pattern,
        prefix=prefix,
        args=args,
    )
    found = await results.get_file(f"{prefix}.tsv")
    if found is None:
        raise FileNotFoundError(f"mmseqs createtsv wrote no {prefix}.tsv")
    return found


async def mmseqs_taxonomy(
    query: Dir,
    target: Dir,
    prefix: str = "sample",
    query_pattern: str = DEFAULT_PATTERN,
    target_pattern: str = DEFAULT_PATTERN,
    args: str = "",
) -> Dir:
    """Taxonomy of the database in ``query`` against the one in ``target``: a Dir of ``<prefix>.*``."""
    return await mmseqs_taxonomy_cmd(
        query=query,
        target=target,
        query_pattern=query_pattern,
        target_pattern=target_pattern,
        prefix=prefix,
        args=args,
    )


async def mmseqs_databases(database: str, args: str = "") -> Dir:
    """Download and build the public database ``database``: a Dir of ``database*``."""
    return await mmseqs_databases_cmd(database=database, args=args)


async def mmseqs_createindex(db: Dir, db_pattern: str = DEFAULT_PATTERN, args: str = "") -> Dir:
    """``db`` with the index of its database (found by ``db_pattern``) added: ``<name>.idx*``."""
    return await mmseqs_createindex_cmd(db=db, db_pattern=db_pattern, args=args)


async def mmseqs_tsv2exprofiledb(database: Dir, args: str = "") -> Dir:
    """``database`` with the expandable profile database of its ``<name>*.tsv`` files added: ``<name>_db*``."""
    return await mmseqs_tsv2exprofiledb_cmd(database=database, args=args)
