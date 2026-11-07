"""duckdb — convert a delimited table to Parquet with DuckDB.

Exposes :func:`duckdb_table2parquet`, which reads a CSV or TSV (optionally
``.gz`` or ``.zst`` compressed) with DuckDB's ``read_csv`` and writes it as
``<prefix>.parquet``. As upstream:

- the delimiter follows the name: ``.csv`` is comma-separated, ``.tsv``
  tab-separated; any other name is an error;
- the first row is the header, and DuckDB infers the column types;
- ``read_options`` (upstream's ``ext.args``) are extra ``read_csv`` options,
  e.g. ``nullstr=['NA','']``, and ``copy_options`` (``ext.args2``) extra
  ``COPY`` options, e.g. ``COMPRESSION 'zstd'``. Both are SQL, passed whole,
  so they may contain spaces.

The table and the output are read and written under fixed names in a work
directory, so neither the input's name nor ``prefix`` ever reaches the SQL. A
``prefix`` containing a newline is refused: Flyte's object store can't hold
such a file name. DuckDB runs on the task's CPUs and memory limit.
"""

import re

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DUCKDB_IMAGE = "community.wave.seqera.io/library/duckdb-cli:1.5.5--9c6d18d9f687a45d"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# The table's name: a .csv or .tsv stem, optionally compressed.
TABLE_NAME = re.compile(r".+\.(csv|tsv)(\.gz|\.zst)?")

# The image puts duckdb on PATH only in its entrypoint, which a shell task
# bypasses, hence the PATH line. The heredoc is unquoted so the limits and
# options expand; expanded values are not expanded again.
duckdb_table2parquet_cmd = shell.create(
    name="duckdb_table2parquet",
    image=DUCKDB_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "table": File,
        "table_name": str,
        "delimiter": str,
        "prefix": str,
        "read_options": str,
        "copy_options": str,
    },
    defaults={"read_options": "", "copy_options": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        export PATH="/opt/wave/.pixi/envs/default/bin:$PATH"
        IN=({{inputs.table}})
        NAME={{inputs.table_name}}
        DELIM={{inputs.delimiter}}
        P={{inputs.prefix}}
        READ={{inputs.read_options}}
        COPY={{inputs.copy_options}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        if [ -n "$READ" ]; then READ=", $READ"; fi
        if [ -n "$COPY" ]; then COPY=", $COPY"; fi
        ln -s "${{IN[0]}}" "table.$NAME"
        duckdb <<SQL
SET threads=$TASK_CPUS;
SET memory_limit='${{TASK_MEMORY_MB}}MB';
SET temp_directory='.';
COPY (SELECT * FROM read_csv('table.$NAME', delim='$DELIM', header=true$READ)) TO 'out.parquet' (FORMAT PARQUET$COPY);
SQL
        mv out.parquet "$OUT/$P.parquet"
    """,
)

env = flyte.TaskEnvironment.from_task("duckdb", duckdb_table2parquet_cmd.as_task())


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def duckdb_table2parquet(
    table: File, prefix: str = "sample", read_options: str = "", copy_options: str = ""
) -> File:
    """``table`` (a ``.csv`` or ``.tsv``, optionally ``.gz``/``.zst``) as ``<prefix>.parquet``."""
    if "\n" in prefix:
        raise ValueError("duckdb_table2parquet: the prefix can't contain a newline (object stores can't hold the name)")
    name = file_name(table)
    match = TABLE_NAME.fullmatch(name)
    if match is None:
        raise ValueError(
            f"DUCKDB_TABLE2PARQUET: cannot determine a delimiter for '{name}' -- "
            "expected a .csv or .tsv suffix, optionally followed by .gz or .zst"
        )
    # The work-directory name keeps the suffix, from which DuckDB picks the compression.
    suffix = name[match.start(1) :]
    delimiter = r"\t" if match.group(1) == "tsv" else ","
    results = await duckdb_table2parquet_cmd(
        table=table,
        table_name=suffix,
        delimiter=delimiter,
        prefix=prefix,
        read_options=read_options,
        copy_options=copy_options,
    )
    found = await results.get_file(f"{prefix}.parquet")
    if found is None:
        raise FileNotFoundError(f"duckdb wrote no {prefix}.parquet")
    return found
