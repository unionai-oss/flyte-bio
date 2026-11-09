"""eggnogmapper — functional annotation of proteins with eggNOG-mapper.

Exposes :func:`eggnogmapper`, ``emapper.py`` on a protein FASTA (may be
gzipped): it finds each protein's seed ortholog in the eggNOG database, then
annotates it from the ortholog's groups (COG category, description, GO terms,
EC, KEGG, CAZy, PFAMs, ...). Outputs keep eggNOG-mapper's names under
``<prefix>``:

- ``<prefix>.emapper.annotations``: the annotation table (always);
- ``<prefix>.emapper.seed_orthologs``: each protein's seed ortholog;
- ``<prefix>.emapper.hits``: the search hits (the last two aren't written by
  every mode, e.g. not by ``no_search``).

``search_mode`` is ``emapper.py -m``. ``db`` is passed with the flag that mode
takes, as upstream: ``diamond`` and ``novel_fams`` a DIAMOND ``.dmnd`` file
(``--dmnd_db``), ``mmseqs`` an MMseqs2 database Dir (``--mmseqs_db``, the
database found by its ``*.dbtype`` file), ``hmmer`` a database
(``--database``), ``no_search`` a ``.seed_orthologs`` table to annotate
(``--annotate_hits_table``) and ``cache`` an annotation cache (``--cache``).
Without ``db`` eggNOG-mapper uses the databases in ``data_dir``.

``data_dir`` is eggNOG-mapper's ``--data_dir``: a Dir holding at least
``eggnog.db`` (the annotation database), plus whatever databases the search
mode reads from it (``eggnog_proteins.dmnd``, ``mmseqs/``, ``pfam/``). As
upstream, eggNOG-mapper runs with ``--cpu`` of the task's CPUs and loads the
annotation database into memory (``--dbmem``) when the task has more than
40000 MiB.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

EGGNOGMAPPER_IMAGE = "quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2"
# Upstream's process_high. The full eggNOG database (tens of GiB) needs the
# task's disk raised, and more than 40000 MiB of memory to load it with --dbmem.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
DB_FLAGS = {
    "diamond": "--dmnd_db",
    "novel_fams": "--dmnd_db",
    "mmseqs": "--mmseqs_db",
    "hmmer": "--database",
    "no_search": "--annotate_hits_table",
    "cache": "--cache",
}

eggnogmapper_cmd = shell.create(
    name="eggnogmapper",
    image=EGGNOGMAPPER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "search_mode": str,
        "db_flag": str,
        "db": File | None,
        "db_dir": Dir | None,
        "data_dir": Dir,
        "prefix": str,
        "args": str,
    },
    defaults={"db_flag": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; DB=({{inputs.db}}); shopt -u nullglob
        DB_DIR={{inputs.db_dir}}
        DATA_DIR={{inputs.data_dir}}
        MODE={{inputs.search_mode}}
        DB_FLAG={{inputs.db_flag}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        if [ "${{NAME##*.}}" = gz ]; then
            gzip -c -d "${{FA[0]}}" > "${{NAME%.gz}}"
            INPUT="${{NAME%.gz}}"
        else
            INPUT="${{FA[0]}}"
        fi
        DB_ARG=()
        if [ -n "$DB_FLAG" ] && [ ${{#DB[@]}} -gt 0 ]; then
            DB_ARG=("$DB_FLAG" "${{DB[0]}}")
        elif [ -n "$DB_FLAG" ] && [ -n "$DB_DIR" ]; then
            # The database's base name: the first *.dbtype file's, as MMseqs2 names them.
            MARKER=$(find -L "$DB_DIR" -maxdepth 1 -name '*.dbtype' | sort | head -n 1)
            if [ -z "$MARKER" ]; then echo "no database in $DB_DIR" >&2; exit 1; fi
            DB_ARG=("$DB_FLAG" "${{MARKER%.dbtype}}")
        fi
        DBMEM=()
        if [ "$TASK_MEMORY_MB" -gt 40000 ]; then DBMEM=(--dbmem); fi
        emapper.py \
            $ARGS \
            --cpu $TASK_CPUS \
            -i "$INPUT" \
            --data_dir "$DATA_DIR" \
            -m "$MODE" \
            "${{DB_ARG[@]}}" \
            "${{DBMEM[@]}}" \
            --output "$P" || exit 1
        mv "$P".emapper.* "$OUT/"
    """,
)

env = flyte.TaskEnvironment.from_task("eggnogmapper", eggnogmapper_cmd.as_task())


@dataclass
class EggnogmapperResult:
    annotations: File  # <prefix>.emapper.annotations
    orthologs: File | None  # <prefix>.emapper.seed_orthologs
    hits: File | None  # <prefix>.emapper.hits
    results: Dir  # everything eggNOG-mapper wrote


async def eggnogmapper(
    fasta: File,
    search_mode: str,
    data_dir: Dir,
    db: File | Dir | None = None,
    prefix: str = "sample",
    args: str = "",
) -> EggnogmapperResult:
    """eggNOG-mapper annotations of the proteins in ``fasta`` (see the module docstring for ``db``)."""
    results = await eggnogmapper_cmd(
        fasta=fasta,
        search_mode=search_mode,
        db_flag=DB_FLAGS.get(search_mode, "") if db is not None else "",
        db=db if isinstance(db, File) else None,
        db_dir=db if isinstance(db, Dir) else None,
        data_dir=data_dir,
        prefix=prefix,
        args=args,
    )
    annotations = await results.get_file(f"{prefix}.emapper.annotations")
    if annotations is None:
        raise FileNotFoundError(f"eggNOG-mapper wrote no {prefix}.emapper.annotations")
    return EggnogmapperResult(
        annotations=annotations,
        orthologs=await results.get_file(f"{prefix}.emapper.seed_orthologs"),
        hits=await results.get_file(f"{prefix}.emapper.hits"),
        results=results,
    )
