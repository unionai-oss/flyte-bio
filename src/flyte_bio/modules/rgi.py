"""rgi — RGI, the Resistance Gene Identifier: antimicrobial resistance genes from the CARD database.

Exposes:

- :func:`rgi_cardannotation`, which prepares a CARD data download (a tar
  archive of ``card.json`` and its companion files, any compression) for RGI:
  ``rgi card_annotation`` writes the reference FASTAs
  ``card_database_v<version>.fasta`` and ``card_database_v<version>_all.fasta``,
  and, as upstream, the processed database is those two FASTAs plus all of the
  download's files. It returns that database as a Dir with its CARD version.
- :func:`rgi_bwt`, which maps one sample's reads (single- or paired-end) to
  CARD's reference alleles (``rgi bwt``, KMA by default). It writes, under
  ``<prefix>.``, the allele mapping JSON (``*.json``), the mapping tables and
  statistics (``*.txt``, among them ``<prefix>.gene_mapping_data.txt``), and
  intermediate files under ``temp/``.
- :func:`rgi_main`, which predicts resistance genes in contigs or a genome
  (``rgi main``): ``<prefix>.json``, ``<prefix>.txt`` and intermediate files
  under ``temp/``.

Both analyses first load the processed database into a working-directory copy
(``rgi load --local``), with the WildCARD variants database too when one is
given (a Dir holding ``wildcard_database_v<version>.fasta``, its ``_all``
twin and ``wildcard/`` with the 61-mer index). As upstream, ``args`` go to
``rgi load`` and ``args2`` to ``rgi bwt`` or ``rgi main``, which run on the
task's CPUs (``--threads``). Intermediate files (``*.xml``, ``*.fsa``, BLAST
indexes, ``*.draft``, ``*.potentialGenes`` and the per-model JSONs) are moved
to ``temp/``, as upstream.

CARD's data is free for academic use only: its README prohibits use by
commercial organisations without a licence from McMaster University.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RGI_IMAGE = "community.wave.seqera.io/library/rgi_kma:e905ecb8305e2609"
# Upstream labels all three process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Sets DB_VERSION from the processed database $CARD, then loads it (and the
# WildCARD database $WILDCARD, when set) into ./localDB, passing $ARGS.
LOAD_DATABASE = r"""
        ALL=("$CARD"/card_database_v*_all.fasta)
        if [ ! -e "${ALL[0]}" ]; then echo "no card_database_v*_all.fasta in the database" >&2; exit 1; fi
        DB_VERSION=$(basename "${ALL[0]}" | sed 's/^card_database_v\([0-9].*[0-9]\).*/\1/')
        WILDCARD_ARGS=()
        if [ -n "$WILDCARD" ]; then
            WILDCARD_ARGS=(
                --wildcard_annotation "$WILDCARD/wildcard_database_v$DB_VERSION.fasta"
                --wildcard_annotation_all_models "$WILDCARD/wildcard_database_v${DB_VERSION}_all.fasta"
                --wildcard_index "$WILDCARD/wildcard/index-for-model-sequences.txt"
                --amr_kmers "$WILDCARD/wildcard/all_amr_61mers.txt"
                --kmer_database "$WILDCARD/wildcard/61_kmer_db.json"
                --kmer_size 61
            )
        fi
        rgi load \
            $ARGS \
            --local \
            --card_json "$CARD/card.json" \
            --debug \
            --card_annotation "$CARD/card_database_v$DB_VERSION.fasta" \
            --card_annotation_all_models "$CARD/card_database_v${DB_VERSION}_all.fasta" \
            "${WILDCARD_ARGS[@]}" || exit 1
"""

# Moves intermediate files to temp/, then the outputs to $OUT.
COLLECT_OUTPUTS = r"""
        mkdir temp
        shopt -s nullglob
        for FILE in *.xml *.fsa *.{nhr,nin,nsq} *.draft *.potentialGenes \
            *{variant,rrna,protein,predictedGenes,overexpression,homolog}.json; do
            mv "$FILE" temp/
        done
        mv ./*.json ./*.txt temp "$OUT"/
"""

rgi_cardannotation_cmd = shell.create(
    name="rgi_cardannotation",
    image=RGI_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"card": File, "args": str},
    defaults={"args": ""},
    outputs={"db": Dir},
    script=r"""
        CARD_ARCHIVE=({inputs.card})
        ARGS={inputs.args}
        OUT={outputs.db}
        W=$(mktemp -d); cd "$W"
        mkdir download
        python -c 'import sys, tarfile; tarfile.open(sys.argv[1]).extractall(sys.argv[2])' \
            "${CARD_ARCHIVE[0]}" download || exit 1
        ENTRIES=(download/*)
        CARD=download
        if [ ${#ENTRIES[@]} -eq 1 ] && [ -d "${ENTRIES[0]}" ]; then CARD="${ENTRIES[0]}"; fi
        rgi card_annotation \
            -i "$CARD/card.json" \
            $ARGS || exit 1
        mv card*.fasta "$OUT"/
        cp -L "$CARD"/* "$OUT"/
    """,
)

rgi_bwt_cmd = shell.create(
    name="rgi_bwt",
    image=RGI_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "card": Dir, "wildcard": Dir | None, "prefix": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        CARD={{inputs.card}}
        WILDCARD={{inputs.wildcard}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {LOAD_DATABASE}
        READ_TWO=()
        if [ ${{#READS[@]}} -gt 1 ]; then READ_TWO=(--read_two "${{READS[1]}}"); fi
        rgi bwt \
            $ARGS2 \
            --local \
            --threads $TASK_CPUS \
            --output_file "$P" \
            --read_one "${{READS[0]}}" \
            "${{READ_TWO[@]}}" || exit 1
        {COLLECT_OUTPUTS}
    """,
)

rgi_main_cmd = shell.create(
    name="rgi_main",
    image=RGI_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "card": Dir, "wildcard": Dir | None, "prefix": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        CARD={{inputs.card}}
        WILDCARD={{inputs.wildcard}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export MPLCONFIGDIR="$W"
        {LOAD_DATABASE}
        ln -s "${{FA[0]}}" .
        rgi main \
            $ARGS2 \
            --threads $TASK_CPUS \
            --output_file "$P" \
            --input_sequence "$(basename "${{FA[0]}}")" || exit 1
        {COLLECT_OUTPUTS}
    """,
)

env = flyte.TaskEnvironment.from_task(
    "rgi", rgi_cardannotation_cmd.as_task(), rgi_bwt_cmd.as_task(), rgi_main_cmd.as_task()
)


@dataclass
class RgiDatabase:
    db: Dir  # the processed CARD database: card_database_v<version>{,_all}.fasta and the download's files
    db_version: str  # CARD's version, e.g. "4.0.2"


@dataclass
class RgiBwtResult:
    json: list[File]  # <prefix>.*.json (the allele mapping)
    tsv: list[File]  # <prefix>.*.txt: mapping tables and statistics
    temp: list[File]  # intermediate files, under temp/
    results: Dir  # everything above


@dataclass
class RgiMainResult:
    json: File  # <prefix>.json
    tsv: File  # <prefix>.txt
    temp: list[File]  # intermediate files, under temp/
    results: Dir  # everything above


async def rgi_cardannotation(card: File, args: str = "") -> RgiDatabase:
    """The CARD database in ``card`` (a data download archive), processed for RGI."""
    db = await rgi_cardannotation_cmd(card=card, args=args)
    return await rgi_database(db)


async def rgi_database(db: Dir) -> RgiDatabase:
    """``db`` with the CARD version its ``card_database_v<version>_all.fasta`` names."""
    names = [f.name for f in await db.list_files()]
    versions = [n.removeprefix("card_database_v").removesuffix("_all.fasta") for n in names if n.endswith("_all.fasta")]
    if len(versions) != 1:
        raise FileNotFoundError(f"expected one card_database_v*_all.fasta in the database, got {names}")
    return RgiDatabase(db=db, db_version=versions[0])


async def split_outputs(results: Dir) -> tuple[list[File], list[File], list[File]]:
    """``results``' top-level JSON and text files, and its temp/ files, each sorted by name."""
    top: list[File] = []
    temp: list[File] = []
    root = results.path.rstrip("/") + "/"
    async for f in results.walk():
        (temp if f.path.removeprefix(root).startswith("temp/") else top).append(f)

    def by_name(files: list[File]) -> list[File]:
        return sorted(files, key=lambda f: f.name)

    json = by_name([f for f in top if f.name.endswith(".json")])
    tsv = by_name([f for f in top if f.name.endswith(".txt")])
    return json, tsv, by_name(temp)


async def rgi_bwt(
    reads: list[File],
    card: Dir,
    wildcard: Dir | None = None,
    prefix: str = "sample",
    args: str = "",
    args2: str = "",
) -> RgiBwtResult:
    """Map one sample's ``reads`` (one or two FASTQs) to the processed CARD database ``card``."""
    results = await rgi_bwt_cmd(reads=reads, card=card, wildcard=wildcard, prefix=prefix, args=args, args2=args2)
    return await rgi_bwt_result(results)


async def rgi_bwt_result(results: Dir) -> RgiBwtResult:
    """Pick ``rgi bwt``'s outputs out of ``results``."""
    json, tsv, temp = await split_outputs(results)
    return RgiBwtResult(json=json, tsv=tsv, temp=temp, results=results)


async def rgi_main(
    fasta: File,
    card: Dir,
    wildcard: Dir | None = None,
    prefix: str = "sample",
    args: str = "",
    args2: str = "",
) -> RgiMainResult:
    """Resistance genes predicted in ``fasta`` (contigs or a genome) with the processed CARD database ``card``."""
    results = await rgi_main_cmd(fasta=fasta, card=card, wildcard=wildcard, prefix=prefix, args=args, args2=args2)
    return await rgi_main_result(results, prefix)


async def rgi_main_result(results: Dir, prefix: str) -> RgiMainResult:
    """Pick ``rgi main``'s outputs out of ``results``."""
    json, tsv, temp = await split_outputs(results)
    json_file = next((f for f in json if f.name == f"{prefix}.json"), None)
    tsv_file = next((f for f in tsv if f.name == f"{prefix}.txt"), None)
    if json_file is None or tsv_file is None:
        raise FileNotFoundError(f"rgi main wrote no {prefix}.json or {prefix}.txt")
    return RgiMainResult(json=json_file, tsv=tsv_file, temp=temp, results=results)
