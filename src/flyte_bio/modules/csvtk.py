"""csvtk — join and concatenate CSV/TSV tables.

Exposes:

- :func:`csvtk_join` — ``csvtk join`` of two or more tables, in list order
  (the first is the left table), on the ``fields`` key(s), written to
  ``<prefix>.tsv`` when the output is tab-separated (``out_delimiter`` a tab,
  or ``-T``/``--out-tabs`` in ``args``) and ``<prefix>.csv`` otherwise.
- :func:`csvtk_concat` — ``csvtk concat`` of tables, in list order, read in
  ``in_format`` and written in ``out_format`` (``tsv``, ``csv``, or a literal
  delimiter) to ``<prefix>.tsv`` for ``tsv`` and ``<prefix>.csv`` otherwise.
  As upstream, the same table may be listed more than once.

Both run on the task's CPUs (``--num-cpus``). The values that upstream quotes
inside ``ext.args`` (the key fields, the delimiters) are their own inputs here,
since ``args`` is split on whitespace without quote removal: write
``fields="ID;ID"`` and ``delimiter="\\t"`` rather than ``--fields 'ID;ID'
-d "\\t"``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CSVTK_IMAGE = "community.wave.seqera.io/library/csvtk:0.37.0--113625988dd3285d"
# Upstream's process_single and process_low; tables are small.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
CONCAT_RESOURCES = flyte.Resources(cpu=2, memory="2Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
csvtk_join_cmd = shell.create(
    name="csvtk_join",
    image=CSVTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "tables": list[File],
        "names": str,
        "fields": str,
        "delimiter": str,
        "out_delimiter": str,
        "filename": str,
        "args": str,
    },
    defaults={"fields": "", "delimiter": "", "out_delimiter": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.tables}})
        NAMES={{inputs.names}}
        FIELDS={{inputs.fields}}
        DELIM={{inputs.delimiter}}
        OUT_DELIM={{inputs.out_delimiter}}
        NAME={{inputs.filename}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        W=$(mktemp -d); cd "$W"
        D=$(dirname "${{IN[0]}}")
        FILES=()
        while IFS= read -r TABLE; do
            ln -s "$D/$TABLE" "$TABLE"
            FILES+=("$TABLE")
        done <<< "$NAMES"
        OPTS=()
        if [ -n "$FIELDS" ]; then OPTS+=(--fields "$FIELDS"); fi
        if [ -n "$DELIM" ]; then OPTS+=(--delimiter "$DELIM"); fi
        if [ -n "$OUT_DELIM" ]; then OPTS+=(--out-delimiter "$OUT_DELIM"); fi
        csvtk join $ARGS "${{OPTS[@]}}" --num-cpus $TASK_CPUS --out-file "$OUT/$NAME" "${{FILES[@]}}" || exit 1
    """,
)

# Each listed table gets its own directory (as upstream's ``inputs/csv*/``),
# so a table listed twice is read twice.
csvtk_concat_cmd = shell.create(
    name="csvtk_concat",
    image=CSVTK_IMAGE,
    resources=CONCAT_RESOURCES,
    inputs={
        "tables": list[File],
        "names": str,
        "delimiter": str,
        "out_delimiter": str,
        "filename": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.tables}})
        NAMES={{inputs.names}}
        DELIM={{inputs.delimiter}}
        OUT_DELIM={{inputs.out_delimiter}}
        NAME={{inputs.filename}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        W=$(mktemp -d); cd "$W"
        D=$(dirname "${{IN[0]}}")
        FILES=()
        N=0
        while IFS= read -r TABLE; do
            N=$((N + 1))
            mkdir -p "inputs/csv$N"
            ln -s "$D/$TABLE" "inputs/csv$N/$TABLE"
            FILES+=("inputs/csv$N/$TABLE")
        done <<< "$NAMES"
        csvtk concat $ARGS --num-cpus $TASK_CPUS --delimiter "$DELIM" --out-delimiter "$OUT_DELIM" \
            --out-file "$OUT/$NAME" "${{FILES[@]}}" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("csvtk", csvtk_join_cmd.as_task(), csvtk_concat_cmd.as_task())


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


def staged(tables: list[File], tool: str) -> tuple[list[File], str]:
    """The distinct tables to stage, and every table's name in list order.

    The same file may be listed more than once; two different files with one
    name can't be staged side by side.
    """
    distinct: dict[str, File] = {}
    for table in tables:
        distinct.setdefault(table.path, table)
    names = [file_name(t) for t in distinct.values()]
    if len(set(names)) != len(names):
        raise ValueError(f"{tool} inputs need distinct names, got {names}")
    return list(distinct.values()), "\n".join(file_name(t) for t in tables)


def format_delimiter(table_format: str) -> str:
    """Upstream's delimiter for a format: a tab for ``tsv``, a comma for ``csv``, else the value itself."""
    return {"tsv": "\t", "csv": ","}.get(table_format, table_format)


async def csvtk_join(
    tables: list[File],
    prefix: str = "sample",
    fields: str = "",
    delimiter: str = "",
    out_delimiter: str = "",
    args: str = "",
) -> File:
    """``csvtk join`` of ``tables`` in list order: ``<prefix>.tsv`` for tab-separated output, else ``<prefix>.csv``.

    ``fields`` is csvtk's ``--fields`` (``ID;ID`` joins on each table's ``ID``);
    ``delimiter`` and ``out_delimiter`` are the input and output delimiters,
    csvtk's defaults when empty.
    """
    if len(tables) < 2:
        raise ValueError("csvtk_join needs at least two tables")
    if len({t.path for t in tables}) != len(tables):
        raise ValueError("csvtk_join needs distinct tables")
    files, names = staged(tables, "csvtk_join")
    words = args.split()
    tabs = out_delimiter == "\t" or "-T" in words or "--out-tabs" in words
    filename = f"{prefix}.{'tsv' if tabs else 'csv'}"
    results = await csvtk_join_cmd(
        tables=files,
        names=names,
        fields=fields,
        delimiter=delimiter,
        out_delimiter=out_delimiter,
        filename=filename,
        args=args,
    )
    return await pick(results, filename, "csvtk join")


async def csvtk_concat(
    tables: list[File], in_format: str, out_format: str, prefix: str = "sample", args: str = ""
) -> File:
    """``csvtk concat`` of ``tables`` in list order: ``<prefix>.tsv`` for ``tsv`` output, else ``<prefix>.csv``."""
    if not tables:
        raise ValueError("csvtk_concat needs at least one table")
    files, names = staged(tables, "csvtk_concat")
    filename = f"{prefix}.{'tsv' if out_format == 'tsv' else 'csv'}"
    results = await csvtk_concat_cmd(
        tables=files,
        names=names,
        delimiter=format_delimiter(in_format),
        out_delimiter=format_delimiter(out_format),
        filename=filename,
        args=args,
    )
    return await pick(results, filename, "csvtk concat")
