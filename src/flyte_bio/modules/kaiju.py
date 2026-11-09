"""kaiju — protein-level taxonomic classification of reads, and its databases.

Exposes:

- :func:`kaiju_mkfmi`, which builds a Kaiju index from a protein FASTA
  (``kaiju-mkbwt``, then ``kaiju-mkfmi``). The database is a Dir holding, as
  upstream's outputs, ``<prefix>.fmi`` with the taxonomy it goes with
  (``nodes.dmp`` and ``names.dmp``, copied in), plus, with
  ``keep_intermediate``, the Burrows-Wheeler transform ``<prefix>.bwt`` and
  suffix array ``<prefix>.sa``. ``-n`` is the task's CPUs; options such as
  the alphabet (``-a``) go in ``args``.
- :func:`kaiju`, which classifies one sample's reads (single- or paired-end)
  against a database Dir (an ``*.fmi`` with ``nodes.dmp`` and ``names.dmp``,
  as :func:`kaiju_mkfmi` builds or a downloaded database extracts to) and
  returns the per-read classifications ``<prefix>.tsv``. ``-z`` is the
  task's CPUs.
- :func:`kaiju2table`, which summarises classifications (one or more
  :func:`kaiju` outputs) at a taxonomic rank (``species``, ``genus``, ...)
  as the table ``<prefix>.txt``, whose first column names each input file.
- :func:`kaiju2krona`, which converts classifications into Krona's text
  format, ``<prefix>.txt``, for ``ktImportText``.

The taxonomy files are found in the database Dir by name (``*nodes.dmp``,
``*names.dmp``), as upstream.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KAIJU_IMAGE = "quay.io/biocontainers/kaiju:1.10.0--h43eeafb_0"
# Upstream labels it process_high; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

kaiju_mkfmi_cmd = shell.create(
    name="kaiju_mkfmi",
    image=KAIJU_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "nodes_dmp": File,
        "names_dmp": File,
        "keep_intermediate": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"keep_intermediate": False, "args": ""},
    outputs={"db": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}}); NODES=({{inputs.nodes_dmp}}); NAMES=({{inputs.names_dmp}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.db}}
        W=$(mktemp -d); cd "$W"
        kaiju-mkbwt \
            $ARGS \
            -n $TASK_CPUS \
            -o "$P" \
            "${{FASTA[0]}}" \
            || exit 1
        kaiju-mkfmi "$P" || exit 1
        cp "${{NODES[0]}}" "$OUT/nodes.dmp"
        cp "${{NAMES[0]}}" "$OUT/names.dmp"
        mv "$P.fmi" "$OUT"/
        if [ {{inputs.keep_intermediate}} = true ]; then mv "$P.bwt" "$P.sa" "$OUT"/; fi
        true
    """,
)

# Upstream labels it process_high; real databases need far more memory.
CLASSIFY_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
# Upstream: process_single.
TABLE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# The database's files, found by name as upstream: its taxonomy tree and names, and its index.
DATABASE_FILES = r"""
        DB_NODES=$(find -L "$DB_DIR" -name "*nodes.dmp" | head -n 1)
        DB_NAMES=$(find -L "$DB_DIR" -name "*names.dmp" | head -n 1)
        DB_INDEX=$(find -L "$DB_DIR" -name "*.fmi" -not -name "._*" | head -n 1)
""".strip()

kaiju_cmd = shell.create(
    name="kaiju",
    image=KAIJU_IMAGE,
    resources=CLASSIFY_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        DB_DIR={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {DATABASE_FILES}
        if [ -z "$DB_NODES" ] || [ -z "$DB_NAMES" ]; then
            echo "the kaiju database is missing one of nodes.dmp, names.dmp" >&2; exit 1
        fi
        if [ -z "$DB_INDEX" ]; then echo "the kaiju database has no .fmi index" >&2; exit 1; fi
        INPUT=(-i "${{R1[0]}}")
        if [ ${{#R2[@]}} -gt 0 ]; then INPUT+=(-j "${{R2[0]}}"); fi
        W=$(mktemp -d); cd "$W"
        kaiju \
            $ARGS \
            -z $TASK_CPUS \
            -t "$DB_NODES" \
            -f "$DB_INDEX" \
            -o "$P.tsv" \
            "${{INPUT[@]}}" \
            || exit 1
        mv "$P.tsv" "$OUT"/
    """,
)

kaiju2table_cmd = shell.create(
    name="kaiju2table",
    image=KAIJU_IMAGE,
    resources=TABLE_RESOURCES,
    inputs={"results": list[File], "db": Dir, "taxon_rank": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"summary": Dir},
    script=rf"""
        IN=({{inputs.results}})
        DB_DIR={{inputs.db}}
        RANK={{inputs.taxon_rank}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.summary}}
        {DATABASE_FILES}
        W=$(mktemp -d); cd "$W"
        # The table's first column names each input file, so they're read by name, as upstream stages them.
        NAMES=()
        for F in "${{IN[@]}}"; do ln -s "$F" "$(basename "$F")"; NAMES+=("$(basename "$F")"); done
        kaiju2table $ARGS \
            -t "$DB_NODES" \
            -n "$DB_NAMES" \
            -r "$RANK" \
            -o "$P.txt" \
            "${{NAMES[@]}}" \
            || exit 1
        mv "$P.txt" "$OUT"/
    """,
)

kaiju2krona_cmd = shell.create(
    name="kaiju2krona",
    image=KAIJU_IMAGE,
    resources=TABLE_RESOURCES,
    inputs={"results": File, "db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"txt": Dir},
    script=rf"""
        IN=({{inputs.results}})
        DB_DIR={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.txt}}
        {DATABASE_FILES}
        W=$(mktemp -d); cd "$W"
        kaiju2krona \
            $ARGS \
            -t "$DB_NODES" \
            -n "$DB_NAMES" \
            -i "${{IN[0]}}" \
            -o "$P.txt" \
            || exit 1
        mv "$P.txt" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "kaiju",
    kaiju_mkfmi_cmd.as_task(),
    kaiju_cmd.as_task(),
    kaiju2table_cmd.as_task(),
    kaiju2krona_cmd.as_task(),
)


async def kaiju_mkfmi(
    fasta: File,
    nodes_dmp: File,
    names_dmp: File,
    keep_intermediate: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """Build a Kaiju index of protein ``fasta``: a Dir of ``<prefix>.fmi``, ``nodes.dmp`` and ``names.dmp``."""
    return await kaiju_mkfmi_cmd(
        fasta=fasta,
        nodes_dmp=nodes_dmp,
        names_dmp=names_dmp,
        keep_intermediate=keep_intermediate,
        prefix=prefix,
        args=args,
    )


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def kaiju(
    reads_1: File,
    db: Dir,
    reads_2: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Classify one sample's reads against the Kaiju database ``db``: ``<prefix>.tsv``."""
    results = await kaiju_cmd(reads_1=reads_1, reads_2=reads_2, db=db, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.tsv", "kaiju")


async def kaiju2table(
    results: File | list[File],
    db: Dir,
    taxon_rank: str,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Kaiju classifications summarised at ``taxon_rank``, using ``db``'s taxonomy: ``<prefix>.txt``."""
    inputs = results if isinstance(results, list) else [results]
    summary = await kaiju2table_cmd(results=inputs, db=db, taxon_rank=taxon_rank, prefix=prefix, args=args)
    return await pick_file(summary, f"{prefix}.txt", "kaiju2table")


async def kaiju2krona(results: File, db: Dir, prefix: str = "sample", args: str = "") -> File:
    """Kaiju classifications in Krona's text format, using ``db``'s taxonomy: ``<prefix>.txt``."""
    txt = await kaiju2krona_cmd(results=results, db=db, prefix=prefix, args=args)
    return await pick_file(txt, f"{prefix}.txt", "kaiju2krona")
