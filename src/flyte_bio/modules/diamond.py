"""diamond — protein sequence clustering and search with DIAMOND.

Exposes:

- :func:`diamond_cluster`, which clusters the proteins of ``db`` (a protein
  FASTA or a DIAMOND ``.dmnd`` database) with ``diamond cluster`` and returns
  its cluster table, ``<prefix>.tsv``: one line per protein, its cluster's
  representative then the protein. As upstream, DIAMOND's memory limit
  (``-M``) is the task's memory and its threads (``-p``) the task's CPUs.
  ``args`` passes the clustering options (``--approx-id``, ``--member-cover``,
  sensitivity modes and so on).
- :func:`diamond_makedb`, a ``.dmnd`` database, ``<prefix>.dmnd``, from a
  protein FASTA (may be gzipped), optionally with taxonomy
  (``--taxonmap``, ``--taxonnodes``, ``--taxonnames``).
- :func:`diamond_blastp`, a protein FASTA (may be gzipped) searched against a
  ``.dmnd`` database: ``<prefix>.<ext>``, where the extension follows the
  output format ``outfmt`` (0 ``blast``, 5 ``xml``, 6 ``txt``, 100 ``daa``,
  101 ``sam``, 102 ``tsv``, 103 ``paf``; anything else falls back to 6), plus
  ``.gz`` when ``args`` has ``--compress 1``.
- :func:`diamond_blastx`, the same for nucleotide queries, translated, with
  the format chosen by its extension (``blast``, ``xml``, ``txt``, ``daa``,
  ``sam``, ``tsv`` or ``paf``; anything else falls back to ``txt``). It also
  returns DIAMOND's log, ``<prefix>.log``.

``blast_columns`` lists the fields of a tabular format (``qseqid qlen ...``).
The searches and makedb run with ``--threads`` of the task's CPUs.

makedb, blastp and blastx use the DIAMOND build their upstream modules pin
(2.2.1), cluster its own (2.2.2), so each build has an environment; :data:`env`
aggregates them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DIAMOND_IMAGE = "community.wave.seqera.io/library/diamond:2.2.2--910ae987965d0d4d"
# Upstream's process_medium asks for 6 CPUs and 36 GB; a proteome of tens of
# thousands of proteins clusters in a few GiB. Raise with
# diamond_cluster_cmd.as_task().override(resources=...) for large databases.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
# The build upstream's makedb, blastp and blastx modules pin.
DIAMOND_SEARCH_IMAGE = "quay.io/biocontainers/diamond:2.2.1--he361c42_0"
# Upstream: makedb process_medium, the searches process_high. Searching a large
# database (such as NR) needs more memory: override the task's resources.
SEARCH_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
BLASTP_EXTENSIONS = {0: "blast", 5: "xml", 6: "txt", 100: "daa", 101: "sam", 102: "tsv", 103: "paf"}
BLASTX_FORMATS = {ext: outfmt for outfmt, ext in BLASTP_EXTENSIONS.items()}

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

diamond_makedb_cmd = shell.create(
    name="diamond_makedb",
    image=DIAMOND_SEARCH_IMAGE,
    resources=SEARCH_RESOURCES,
    inputs={
        "fasta": File,
        "taxonmap": File | None,
        "taxonnodes": File | None,
        "taxonnames": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob
        TAXONMAP=({{inputs.taxonmap}}); TAXONNODES=({{inputs.taxonnodes}}); TAXONNAMES=({{inputs.taxonnames}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        if [ "${{NAME##*.}}" = gz ]; then
            gzip -c -d "${{FA[0]}}" > "${{NAME%.gz}}"
            IN="${{NAME%.gz}}"
        else
            IN="${{FA[0]}}"
        fi
        TAXONOMY=()
        if [ ${{#TAXONMAP[@]}} -gt 0 ]; then TAXONOMY+=(--taxonmap "${{TAXONMAP[0]}}"); fi
        if [ ${{#TAXONNODES[@]}} -gt 0 ]; then TAXONOMY+=(--taxonnodes "${{TAXONNODES[0]}}"); fi
        if [ ${{#TAXONNAMES[@]}} -gt 0 ]; then TAXONOMY+=(--taxonnames "${{TAXONNAMES[0]}}"); fi
        diamond \
            makedb \
            --threads $TASK_CPUS \
            --in "$IN" \
            -d "$OUT/$P" \
            $ARGS \
            "${{TAXONOMY[@]}}" || exit 1
    """,
)

# As upstream, blastp reads the query as it is (DIAMOND reads gzip).
diamond_blastp_cmd = shell.create(
    name="diamond_blastp",
    image=DIAMOND_SEARCH_IMAGE,
    resources=SEARCH_RESOURCES,
    inputs={
        "fasta": File,
        "db": File,
        "outfmt": int,
        "blast_columns": str,
        "out_name": str,
        "args": str,
    },
    defaults={"blast_columns": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        DB=({{inputs.db}})
        COLUMNS={{inputs.blast_columns}}
        NAME={{inputs.out_name}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        diamond \
            blastp \
            --threads $TASK_CPUS \
            --db "${{DB[0]}}" \
            --query "${{FA[0]}}" \
            --outfmt {{inputs.outfmt}} $COLUMNS \
            $ARGS \
            --out "$OUT/$NAME" || exit 1
    """,
)

diamond_blastx_cmd = shell.create(
    name="diamond_blastx",
    image=DIAMOND_SEARCH_IMAGE,
    resources=SEARCH_RESOURCES,
    inputs={
        "fasta": File,
        "db": File,
        "outfmt": int,
        "blast_columns": str,
        "out_name": str,
        "prefix": str,
        "args": str,
    },
    defaults={"blast_columns": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        DB=({{inputs.db}})
        COLUMNS={{inputs.blast_columns}}
        NAME={{inputs.out_name}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        QUERY=$(basename "${{FA[0]}}")
        if [ "${{QUERY##*.}}" = gz ]; then
            gzip -c -d "${{FA[0]}}" > "${{QUERY%.gz}}"
            QUERY="${{QUERY%.gz}}"
        else
            QUERY="${{FA[0]}}"
        fi
        diamond \
            blastx \
            --threads $TASK_CPUS \
            --db "${{DB[0]}}" \
            --query "$QUERY" \
            --outfmt {{inputs.outfmt}} $COLUMNS \
            $ARGS \
            --out "$OUT/$NAME" \
            --log || exit 1
        mv diamond.log "$OUT/$P.log"
    """,
)

diamond_cluster_env = flyte.TaskEnvironment.from_task("diamond_cluster", diamond_cluster_cmd.as_task())
diamond_search_env = flyte.TaskEnvironment.from_task(
    "diamond_search", diamond_makedb_cmd.as_task(), diamond_blastp_cmd.as_task(), diamond_blastx_cmd.as_task()
)
env = flyte.TaskEnvironment(name="diamond", depends_on=[diamond_cluster_env, diamond_search_env])


async def diamond_cluster(db: File, prefix: str = "sample", args: str = "") -> File:
    """``diamond cluster`` of the proteins in ``db`` (FASTA or ``.dmnd``): ``<prefix>.tsv``."""
    results = await diamond_cluster_cmd(db=db, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.tsv")
    if found is None:
        raise FileNotFoundError(f"diamond cluster wrote no {prefix}.tsv")
    return found


async def picked(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"diamond {tool} wrote no {name}")
    return found


async def diamond_makedb(
    fasta: File,
    prefix: str = "sample",
    taxonmap: File | None = None,
    taxonnodes: File | None = None,
    taxonnames: File | None = None,
    args: str = "",
) -> File:
    """A DIAMOND database of the proteins in ``fasta`` (may be gzipped): ``<prefix>.dmnd``."""
    results = await diamond_makedb_cmd(
        fasta=fasta, taxonmap=taxonmap, taxonnodes=taxonnodes, taxonnames=taxonnames, prefix=prefix, args=args
    )
    return await picked(results, f"{prefix}.dmnd", "makedb")


def compressed(args: str) -> bool:
    """Whether ``args`` asks DIAMOND to gzip its output (``--compress 1``), as upstream checks."""
    words = args.split()
    return any(w == "--compress" and words[i + 1 : i + 2] == ["1"] for i, w in enumerate(words))


async def diamond_blastp(
    fasta: File, db: File, prefix: str = "sample", outfmt: int = 6, blast_columns: str = "", args: str = ""
) -> File:
    """``diamond blastp`` of ``fasta`` against ``db`` (``.dmnd``): ``<prefix>.<ext>`` (see the module docstring)."""
    if outfmt not in BLASTP_EXTENSIONS:
        print(f"Unknown output format {outfmt}: selecting DIAMOND default of tabular BLAST output (txt)")
        outfmt = 6
    name = f"{prefix}.{BLASTP_EXTENSIONS[outfmt]}" + (".gz" if compressed(args) else "")
    results = await diamond_blastp_cmd(
        fasta=fasta, db=db, outfmt=outfmt, blast_columns=blast_columns, out_name=name, args=args
    )
    return await picked(results, name, "blastp")


@dataclass
class DiamondBlastxResult:
    hits: File  # <prefix>.<ext>
    log: File  # <prefix>.log


async def diamond_blastx(
    fasta: File, db: File, prefix: str = "sample", out_ext: str = "txt", blast_columns: str = "", args: str = ""
) -> DiamondBlastxResult:
    """``diamond blastx`` of the nucleotide ``fasta`` (may be gzipped) against ``db``: ``<prefix>.<out_ext>``."""
    if out_ext not in BLASTX_FORMATS:
        print(f"Unknown output format {out_ext}: selecting DIAMOND default of tabular BLAST output (txt)")
        out_ext = "txt"
    name = f"{prefix}.{out_ext}"
    results = await diamond_blastx_cmd(
        fasta=fasta,
        db=db,
        outfmt=BLASTX_FORMATS[out_ext],
        blast_columns=blast_columns,
        out_name=name,
        prefix=prefix,
        args=args,
    )
    return DiamondBlastxResult(
        hits=await picked(results, name, "blastx"), log=await picked(results, f"{prefix}.log", "blastx")
    )
