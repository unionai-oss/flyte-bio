"""catpack — CAT_pack (CAT and BAT), taxonomic classification of contigs and metagenome-assembled genomes.

Exposes:

- :func:`catpack_download` — download a CAT_pack source database (``nr`` or
  ``gtdb``): the protein FASTA ``*.<db>.gz``, ``*.names.dmp``,
  ``*.nodes.dmp``, ``*accession2taxid*.gz`` and the download log.
- :func:`catpack_prepare` — build a CAT_pack database from a protein FASTA,
  NCBI-style ``names.dmp``/``nodes.dmp`` and an accession-to-taxid table: a
  database Dir (the DIAMOND database ``*.dmnd``, ``*.fastaid2LCAtaxid``,
  ``*.taxids_with_multiple_offspring``) and a taxonomy Dir (``names.dmp``,
  ``nodes.dmp``).
- :func:`catpack_contigs` — classify contigs (CAT):
  ``<prefix>.ORF2LCA.txt``, ``<prefix>.contig2classification.txt``,
  ``<prefix>.log`` and, unless premade, ``<prefix>.alignment.diamond``,
  ``<prefix>.predicted_proteins.faa`` and ``<prefix>.predicted_proteins.gff``.
- :func:`catpack_bins` — classify bins (BAT): ``<prefix>.ORF2LCA.txt``,
  ``<prefix>.bin2classification.txt``, ``<prefix>.log`` and, unless premade,
  ``<prefix>.concatenated.alignment.diamond``,
  ``<prefix>.concatenated.predicted_proteins.faa`` and
  ``<prefix>.concatenated.predicted_proteins.gff``.
- :func:`catpack_addnames` — add taxon names to a classification:
  ``<prefix>.txt``.
- :func:`catpack_summarise` — summarise a named classification (with the
  contigs, for contig classifications): ``<prefix>.txt``.

As upstream:

- ``prepare``, ``contigs`` and ``bins`` get the task's CPUs as threads
  (``-n``); ``args`` come last.
- CAT_pack doesn't read gzipped input, so gzipped contigs are decompressed
  first; for bins, a ``bin_suffix`` ending in ``.gz`` decompresses every bin
  and passes the suffix without ``.gz``.
- ``add_names`` runs in upstream's CAT_pack 6.0.1 image, the other tools in
  its 6.0 image.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CATPACK_IMAGE = "community.wave.seqera.io/library/cat_gzip:0ab95a62b35744c9"
CATPACK_ADDNAMES_IMAGE = "community.wave.seqera.io/library/cat:6.0.1--45eb3de71f06f426"

# Upstream: process_single + process_long. The nr database is well over 100 GB; raise the disk with
# catpack_download_cmd.override(resources=...) if it doesn't fit.
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="200Gi")
# Upstream labels prepare process_medium + process_high_memory; these fit test-sized databases. A full nr
# database needs far more memory and disk.
PREPARE_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
# Upstream: process_medium. The 4 CPUs match upstream's test CPUs.
CLASSIFY_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
# Upstream: process_single.
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

catpack_download_cmd = shell.create(
    name="catpack_download",
    image=CATPACK_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"db": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        DB={inputs.db}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        CAT_pack \
            download \
            $ARGS \
            --db "$DB" \
            -o "$P"/ || exit 1
        mv "$P"/* "$OUT"/
    """,
)

catpack_prepare_cmd = shell.create(
    name="catpack_prepare",
    image=CATPACK_IMAGE,
    resources=PREPARE_RESOURCES,
    inputs={"db_fasta": File, "names": File, "nodes": File, "acc2tax": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"db": Dir, "taxonomy": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.db_fasta}}); NAMES=({{inputs.names}}); NODES=({{inputs.nodes}}); ACC2TAX=({{inputs.acc2tax}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for F in "${{FA[0]}}" "${{NAMES[0]}}" "${{NODES[0]}}" "${{ACC2TAX[0]}}"; do ln -s "$F" "$(basename "$F")"; done
        CAT_pack prepare \
            -n $TASK_CPUS \
            --db_fasta "$(basename "${{FA[0]}}")" \
            --names "$(basename "${{NAMES[0]}}")" \
            --nodes "$(basename "${{NODES[0]}}")" \
            --acc2tax "$(basename "${{ACC2TAX[0]}}")" \
            --db_dir "$P"/ \
            $ARGS || exit 1
        mv "$P"/db/* {{outputs.db}}/
        mv "$P"/tax/* {{outputs.taxonomy}}/
    """,
)

catpack_contigs_cmd = shell.create(
    name="catpack_contigs",
    image=CATPACK_IMAGE,
    resources=CLASSIFY_RESOURCES,
    inputs={
        "contigs": list[File],
        "database": Dir,
        "taxonomy": Dir,
        "proteins": File | None,
        "diamond_table": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        CONTIGS=({{inputs.contigs}})
        shopt -s nullglob; PROTEINS=({{inputs.proteins}}); TABLE=({{inputs.diamond_table}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # CAT_pack does not support gzipped input, so decompress any gzipped contigs.
        CONTIGS_FASTAS=()
        for F in "${{CONTIGS[@]}}"; do
            NAME=$(basename "$F")
            if [[ "$NAME" == *.gz ]]; then
                gunzip -c "$F" > "${{NAME%.gz}}"
                CONTIGS_FASTAS+=("${{NAME%.gz}}")
            else
                ln -s "$F" "$NAME"
                CONTIGS_FASTAS+=("$NAME")
            fi
        done
        PREMADE=()
        if [ ${{#PROTEINS[@]}} -gt 0 ]; then PREMADE+=(--proteins_fasta "${{PROTEINS[0]}}"); fi
        if [ ${{#TABLE[@]}} -gt 0 ]; then PREMADE+=(--diamond_alignment "${{TABLE[0]}}"); fi

        CAT_pack contigs \
            --nproc $TASK_CPUS \
            --contigs_fasta "${{CONTIGS_FASTAS[@]}}" \
            --database_folder {{inputs.database}} \
            --taxonomy_folder {{inputs.taxonomy}} \
            --out_prefix "$P" \
            "${{PREMADE[@]}}" \
            $ARGS || exit 1
        shopt -s nullglob
        mv "$P".ORF2LCA.txt "$P".contig2classification.txt "$P".log "$P".*diamond \
            "$P".*predicted_proteins.faa "$P".*gff "$OUT"/
    """,
)

catpack_bins_cmd = shell.create(
    name="catpack_bins",
    image=CATPACK_IMAGE,
    resources=CLASSIFY_RESOURCES,
    inputs={
        "bins": list[File],
        "database": Dir,
        "taxonomy": Dir,
        "proteins": File | None,
        "diamond_table": File | None,
        "bin_suffix": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BINS=({{inputs.bins}})
        shopt -s nullglob; PROTEINS=({{inputs.proteins}}); TABLE=({{inputs.diamond_table}}); shopt -u nullglob
        SUFFIX={{inputs.bin_suffix}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir bins
        for F in "${{BINS[@]}}"; do ln -s "$F" bins/"$(basename "$F")"; done

        # CAT_pack does not support gzipped input. If the bin suffix ends in .gz, decompress the bins and
        # strip .gz off the suffix passed to CAT_pack.
        BINS_DIR=bins
        CAT_SUFFIX=$SUFFIX
        if [[ "$SUFFIX" == *.gz ]]; then
            CAT_SUFFIX=${{SUFFIX%.gz}}
            BINS_DIR=input_bins
            mkdir input_bins
            for F in bins/*"$SUFFIX"; do
                gunzip -c "$F" > "input_bins/$(basename "$F" .gz)"
            done
        fi
        PREMADE=()
        if [ ${{#PROTEINS[@]}} -gt 0 ]; then PREMADE+=(-p "${{PROTEINS[0]}}"); fi
        if [ ${{#TABLE[@]}} -gt 0 ]; then PREMADE+=(-d "${{TABLE[0]}}"); fi

        CAT_pack bins \
            -n $TASK_CPUS \
            -b "$BINS_DIR"/ \
            -d {{inputs.database}} \
            -t {{inputs.taxonomy}} \
            -s "$CAT_SUFFIX" \
            "${{PREMADE[@]}}" \
            -o "$P" \
            $ARGS || exit 1
        shopt -s nullglob
        mv "$P".ORF2LCA.txt "$P".bin2classification.txt "$P".log "$P".*diamond \
            "$P".*predicted_proteins.faa "$P".*gff "$OUT"/
    """,
)

catpack_addnames_cmd = shell.create(
    name="catpack_addnames",
    image=CATPACK_ADDNAMES_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"classification": File, "taxonomy": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.classification})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        ln -s "${IN[0]}" "$(basename "${IN[0]}")"
        CAT_pack add_names \
            $ARGS \
            -i "$(basename "${IN[0]}")" \
            -t {inputs.taxonomy} \
            -o "$P".txt || exit 1
        mv "$P".txt {outputs.results}/
    """,
)

catpack_summarise_cmd = shell.create(
    name="catpack_summarise",
    image=CATPACK_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"classification": File, "contigs": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.classification})
        shopt -s nullglob; CONTIGS=({inputs.contigs}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        ln -s "${IN[0]}" "$(basename "${IN[0]}")"
        INSERT_CONTIGS=()
        if [ ${#CONTIGS[@]} -gt 0 ]; then
            ln -s "${CONTIGS[0]}" "$(basename "${CONTIGS[0]}")"
            INSERT_CONTIGS=(-c "$(basename "${CONTIGS[0]}")")
        fi
        CAT_pack summarise \
            $ARGS \
            -i "$(basename "${IN[0]}")" \
            "${INSERT_CONTIGS[@]}" \
            -o "$P".txt || exit 1
        mv "$P".txt {outputs.results}/
    """,
)

# One image per environment: add_names runs in the 6.0.1 image, the other tools in the 6.0 one.
catpack_env = flyte.TaskEnvironment.from_task(
    "catpack",
    catpack_download_cmd.as_task(),
    catpack_prepare_cmd.as_task(),
    catpack_contigs_cmd.as_task(),
    catpack_bins_cmd.as_task(),
    catpack_summarise_cmd.as_task(),
)
catpack_addnames_env = flyte.TaskEnvironment.from_task("catpack_addnames", catpack_addnames_cmd.as_task())
env = flyte.TaskEnvironment(name="catpack_tools", depends_on=[catpack_env, catpack_addnames_env])


def file_name(file: File) -> str:
    return file.name or file.path.rsplit("/", maxsplit=1)[-1]


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"CAT_pack {tool} wrote no {name}")
    return found


async def pick_suffix(results: Dir, suffix: str, tool: str) -> File:
    found = [f for f in await results.list_files() if file_name(f).endswith(suffix)]
    if len(found) != 1:
        raise FileNotFoundError(f"CAT_pack {tool} wrote {len(found)} *{suffix} files, expected one")
    return found[0]


@dataclass
class CatpackDownloadResult:
    fasta: File  # the protein FASTA, *.<db>.gz
    names: File
    nodes: File
    acc2tax: File
    log: File
    results: Dir  # everything above


async def catpack_download(db: str, prefix: str = "db", args: str = "") -> CatpackDownloadResult:
    """Download the CAT_pack source database ``db`` (``nr`` or ``gtdb``)."""
    results = await catpack_download_cmd(db=db, prefix=prefix, args=args)
    acc2tax = [f for f in await results.list_files() if "accession2taxid" in file_name(f)]
    if len(acc2tax) != 1:
        raise FileNotFoundError(f"CAT_pack download wrote {len(acc2tax)} accession2taxid files, expected one")
    return CatpackDownloadResult(
        fasta=await pick_suffix(results, f".{db}.gz", "download"),
        names=await pick_suffix(results, ".names.dmp", "download"),
        nodes=await pick_suffix(results, ".nodes.dmp", "download"),
        acc2tax=acc2tax[0],
        log=await pick_suffix(results, ".log", "download"),
        results=results,
    )


@dataclass
class CatpackDatabase:
    db: Dir  # the DIAMOND database and CAT_pack's LCA tables
    taxonomy: Dir  # names.dmp and nodes.dmp


async def catpack_prepare(
    db_fasta: File, names: File, nodes: File, acc2tax: File, prefix: str = "db", args: str = ""
) -> CatpackDatabase:
    """Build a CAT_pack database and taxonomy from a protein FASTA and NCBI-style taxonomy files."""
    db, taxonomy = await catpack_prepare_cmd(
        db_fasta=db_fasta, names=names, nodes=nodes, acc2tax=acc2tax, prefix=prefix, args=args
    )
    return CatpackDatabase(db=db, taxonomy=taxonomy)


@dataclass
class CatpackClassifyResult:
    orf2lca: File
    classification: File  # <prefix>.contig2classification.txt or <prefix>.bin2classification.txt
    log: File
    diamond: File | None  # the DIAMOND alignment, unless premade
    faa: File | None  # the predicted proteins, unless premade
    gff: File | None  # the predicted genes, unless premade
    results: Dir  # everything above


async def read_classify_results(results: Dir, prefix: str, kind: str, tool: str) -> CatpackClassifyResult:
    files = await results.list_files()

    def optional(suffix: str) -> File | None:
        found = [f for f in files if file_name(f).endswith(suffix)]
        return found[0] if found else None

    return CatpackClassifyResult(
        orf2lca=await pick(results, f"{prefix}.ORF2LCA.txt", tool),
        classification=await pick(results, f"{prefix}.{kind}2classification.txt", tool),
        log=await pick(results, f"{prefix}.log", tool),
        diamond=optional(".diamond"),
        faa=optional(".predicted_proteins.faa"),
        gff=optional(".gff"),
        results=results,
    )


async def catpack_contigs(
    contigs: File | list[File],
    db: CatpackDatabase,
    proteins: File | None = None,
    diamond_table: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> CatpackClassifyResult:
    """Classify contigs (plain or gzipped FASTA), optionally from premade proteins or DIAMOND alignment."""
    results = await catpack_contigs_cmd(
        contigs=contigs if isinstance(contigs, list) else [contigs],
        database=db.db,
        taxonomy=db.taxonomy,
        proteins=proteins,
        diamond_table=diamond_table,
        prefix=prefix,
        args=args,
    )
    return await read_classify_results(results, prefix, "contig", "contigs")


async def catpack_bins(
    bins: File | list[File],
    db: CatpackDatabase,
    bin_suffix: str,
    proteins: File | None = None,
    diamond_table: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> CatpackClassifyResult:
    """Classify bins: each FASTA ending in ``bin_suffix`` (gzipped if it ends in ``.gz``) is one bin."""
    results = await catpack_bins_cmd(
        bins=bins if isinstance(bins, list) else [bins],
        database=db.db,
        taxonomy=db.taxonomy,
        proteins=proteins,
        diamond_table=diamond_table,
        bin_suffix=bin_suffix,
        prefix=prefix,
        args=args,
    )
    return await read_classify_results(results, prefix, "bin", "bins")


async def catpack_addnames(classification: File, taxonomy: Dir, prefix: str = "sample", args: str = "") -> File:
    """``classification`` with taxon names added: ``<prefix>.txt``."""
    if file_name(classification) == f"{prefix}.txt":
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await catpack_addnames_cmd(classification=classification, taxonomy=taxonomy, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.txt", "add_names")


async def catpack_summarise(
    classification: File, contigs: File | None = None, prefix: str = "sample", args: str = ""
) -> File:
    """Summarise a named classification: ``<prefix>.txt``. Contig classifications also need the contigs."""
    if file_name(classification) == f"{prefix}.txt":
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await catpack_summarise_cmd(classification=classification, contigs=contigs, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.txt", "summarise")
