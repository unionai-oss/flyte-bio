"""blast — BLAST databases and nucleotide search.

Exposes:

- :func:`blast_makeblastdb` — a BLAST database from a FASTA (may be gzipped):
  a Dir of ``<fasta name>.*`` (``.nhr``, ``.nin``, ``.nsq``, ... for
  nucleotides; ``.phr``, ``.pin``, ``.psq``, ... for proteins). As with
  makeblastdb itself the default database type is protein: pass
  ``-dbtype nucl`` in ``args`` for nucleotides. With ``taxid_map``, sequences
  get taxonomy IDs (``-taxid_map``).
- :func:`blast_blastn` — search a nucleotide FASTA (may be gzipped) against a
  nucleotide database Dir (a multi-volume database's ``.nal`` alias if it has
  one, else its ``.nin``): ``<prefix>.txt``. Hits can be limited to taxonomy
  IDs, from a file (``taxidlist``) or a comma-separated list (``taxids``), or,
  with ``negative_tax``, exclude them. The output format, which has spaces
  (``6 qseqid sseqid pident``), goes in ``outfmt`` rather than ``args``.

The database records its creation date (``.nin``, ``.njs``), so only its
sequence, header and taxonomy files are byte-reproducible; the default report
records the database's date too.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BLAST_IMAGE = "community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759"
# Upstream labels both process_medium; searching large databases needs more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

blast_makeblastdb_cmd = shell.create(
    name="blast_makeblastdb",
    image=BLAST_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "taxid_map": File | None, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        shopt -s nullglob; TAXID_MAP=({inputs.taxid_map}); shopt -u nullglob
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        NAME=$(basename "${FA[0]}")
        if [ "${NAME##*.}" = gz ]; then
            NAME=${NAME%.gz}
            gzip -c -d "${FA[0]}" > "$NAME"
        else
            ln -s "${FA[0]}" "$NAME"
        fi
        EXTRA=()
        if [ ${#TAXID_MAP[@]} -gt 0 ]; then EXTRA+=(-taxid_map "${TAXID_MAP[0]}"); fi
        # Built under a relative path, as upstream, in case the database records it.
        mkdir db
        makeblastdb -in "$NAME" -out "db/$NAME" $ARGS "${EXTRA[@]}" && mv db/* "$OUT"/
    """,
)

blast_blastn_cmd = shell.create(
    name="blast_blastn",
    image=BLAST_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "db": Dir,
        "taxidlist": File | None,
        "taxids": str,
        "negative_tax": bool,
        "outfmt": str,
        "prefix": str,
        "args": str,
    },
    defaults={"taxids": "", "negative_tax": False, "outfmt": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; TAXIDLIST=({{inputs.taxidlist}}); shopt -u nullglob
        DB_DIR={{inputs.db}}
        TAXIDS={{inputs.taxids}}
        OUTFMT={{inputs.outfmt}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        NAME=$(basename "${{FA[0]}}")
        if [ "${{NAME##*.}}" = gz ]; then
            NAME=${{NAME%.gz}}
            gzip -c -d "${{FA[0]}}" > "$NAME"
        else
            ln -s "${{FA[0]}}" "$NAME"
        fi

        export BLASTDB=$DB_DIR
        DB=$(find -L "$DB_DIR" -name "*.nal" | sed 's/\.nal$//' | head -n 1)
        if [ -z "$DB" ]; then
            DB=$(find -L "$DB_DIR" -name "*.nin" | sed 's/\.nin$//' | head -n 1)
        fi
        if [ -z "$DB" ]; then
            echo "no .nal or .nin database in the db Dir" >&2
            exit 1
        fi
        echo "Using $DB"

        NEGATIVE=""
        if [ {{inputs.negative_tax}} = true ]; then NEGATIVE=negative_; fi
        EXTRA=()
        if [ ${{#TAXIDLIST[@]}} -gt 0 ]; then EXTRA+=(-${{NEGATIVE}}taxidlist "${{TAXIDLIST[0]}}"); fi
        if [ -n "$TAXIDS" ]; then EXTRA+=(-${{NEGATIVE}}taxids "$TAXIDS"); fi
        if [ -n "$OUTFMT" ]; then EXTRA+=(-outfmt "$OUTFMT"); fi

        blastn \
            -num_threads $TASK_CPUS \
            -db "$DB" \
            -query "$NAME" \
            "${{EXTRA[@]}}" \
            $ARGS \
            -out "$OUT/$P.txt"
    """,
)

env = flyte.TaskEnvironment.from_task("blast", blast_makeblastdb_cmd.as_task(), blast_blastn_cmd.as_task())


async def blast_makeblastdb(fasta: File, taxid_map: File | None = None, args: str = "") -> Dir:
    """A BLAST database of ``fasta``: a Dir of ``<fasta name>.*`` (pass ``-dbtype nucl`` for nucleotides)."""
    return await blast_makeblastdb_cmd(fasta=fasta, taxid_map=taxid_map, args=args)


async def blast_blastn(
    fasta: File,
    db: Dir,
    taxidlist: File | None = None,
    taxids: str = "",
    negative_tax: bool = False,
    outfmt: str = "",
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Search ``fasta`` against the nucleotide database ``db``: ``<prefix>.txt``."""
    if taxidlist is not None and taxids:
        raise ValueError("taxidlist and taxids can not be used at the same time, choose only one")
    results = await blast_blastn_cmd(
        fasta=fasta,
        db=db,
        taxidlist=taxidlist,
        taxids=taxids,
        negative_tax=negative_tax,
        outfmt=outfmt,
        prefix=prefix,
        args=args,
    )
    found = await results.get_file(f"{prefix}.txt")
    if found is None:
        raise FileNotFoundError(f"blastn wrote no {prefix}.txt")
    return found
