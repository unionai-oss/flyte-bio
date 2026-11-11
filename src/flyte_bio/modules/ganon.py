"""ganon — taxonomic classification databases (interleaved Bloom filters).

Exposes :func:`ganon_buildcustom`, which builds a ganon database from
reference FASTA files (``ganon build-custom``): either the files themselves
(``--input``), or a TSV that lists them with their targets and taxa
(``--input-file``; it names the files by their basenames). Outputs keep
upstream's names:

- the filter, ``<prefix>.hibf`` (or ``<prefix>.ibf`` with
  ``--filter-type ibf``);
- the taxonomy, ``<prefix>.tax`` (unless ``--taxonomy skip``);
- with ``--write-info-file``, ``<prefix>.info.tsv``.

Without ``taxonomy_files`` ganon downloads the NCBI taxonomy itself, so the
task needs internet access. NCBI taxonomy files are passed in the order ganon
reads them, ``nodes``, ``names`` then ``merged``, whatever order they're given
in. ``--threads`` is the task's CPUs; options such as
``--input-target`` go in ``args``.

The database's filter is found by its ``*.hibf`` / ``*.ibf`` file (its
``.tax`` beside it) by the other tools, which take the database Dir:

- :func:`ganon_classify` classifies one sample's reads (single- or
  paired-end; ``ganon classify``). Outputs keep upstream's names: the report
  ``<prefix>.rep``, the default taxonomic report ``<prefix>.tre``, the log
  ``<prefix>.log`` and, with ``--output-one`` / ``--output-all`` /
  ``--output-unclassified`` in ``args``, the per-read files
  ``<prefix>.{one,all,unc}``. ``--threads`` is the task's CPUs.
- :func:`ganon_report` writes a taxonomic report from a ``.rep``
  (``ganon report``): ``<prefix>.tre``; options such as ``--report-type``
  and ``--output-format`` go in ``args``.
- :func:`ganon_table` combines several samples' ``.tre`` reports into one
  table (``ganon table``), ``<prefix>.txt``; the samples are named after the
  files.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

GANON_IMAGE = "quay.io/biocontainers/ganon:2.1.0--py310hab1bfa5_1"
# Upstream labels it process_high; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

ganon_buildcustom_cmd = shell.create(
    name="ganon_buildcustom",
    image=GANON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": list[File],
        "input_tsv": File | None,
        "taxonomy_files": list[File],
        "genome_size_files": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}})
        shopt -s nullglob
        TSV=({{inputs.input_tsv}}); TAXONOMY=({{inputs.taxonomy_files}}); SIZES=({{inputs.genome_size_files}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${{FASTA[@]}}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        if [ ${{#TSV[@]}} -gt 0 ]; then
            ln -s "${{TSV[0]}}" "$(basename "${{TSV[0]}}")"
            INPUT=(--input-file "$(basename "${{TSV[0]}}")")
        else
            INPUT=(--input "${{NAMES[@]}}")
        fi
        # NCBI taxonomy files go in ganon's order: nodes, names, merged, then any others.
        ORDERED=()
        for PATTERN in '*nodes*' '*names*' '*merged*'; do
            for F in "${{TAXONOMY[@]}}"; do
                case "$(basename "$F")" in $PATTERN) ORDERED+=("$F") ;; esac
            done
        done
        for F in "${{TAXONOMY[@]}}"; do
            case "$(basename "$F")" in *nodes*|*names*|*merged*) ;; *) ORDERED+=("$F") ;; esac
        done
        TAXONOMY_ARGS=()
        if [ ${{#ORDERED[@]}} -gt 0 ]; then TAXONOMY_ARGS=(--taxonomy-files "${{ORDERED[@]}}"); fi
        SIZE_ARGS=()
        if [ ${{#SIZES[@]}} -gt 0 ]; then SIZE_ARGS=(--genome-size-files "${{SIZES[@]}}"); fi
        ganon \
            build-custom \
            --threads $TASK_CPUS \
            "${{INPUT[@]}}" \
            --db-prefix "$P" \
            "${{TAXONOMY_ARGS[@]}}" \
            "${{SIZE_ARGS[@]}}" \
            $ARGS \
            > ganon.log 2>&1 \
            || {{ tail -n 20 ganon.log >&2; exit 1; }}
        cat ganon.log >&2
        for F in "$P".hibf "$P".ibf "$P".tax "$P".info.tsv; do
            if [ -e "$F" ]; then mv "$F" "$OUT"/; fi
        done
    """,
)

# Upstream labels classify process_high, report and table process_single.
CLASSIFY_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# The database's prefix: its filter's path without .hibf / .ibf.
FIND_DB_PREFIX = r"""
        FILTER=$(find -L "$DB" -name '*.*ibf' | head -n 1)
        if [ -z "$FILTER" ]; then echo "no ganon filter (*.hibf, *.ibf) in the database" >&2; exit 1; fi
        DB_PREFIX=$(echo "$FILTER" | sed 's/\.h\?ibf$//')
""".strip()

ganon_classify_cmd = shell.create(
    name="ganon_classify",
    image=GANON_IMAGE,
    resources=CLASSIFY_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        DB={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {FIND_DB_PREFIX}
        W=$(mktemp -d); cd "$W"
        if [ ${{#R2[@]}} -gt 0 ]; then
            INPUT=(--paired-reads "${{R1[0]}}" "${{R2[0]}}")
        else
            INPUT=(--single-reads "${{R1[0]}}")
        fi
        set -o pipefail
        ganon \
            classify \
            --db-prefix "$DB_PREFIX" \
            $ARGS \
            --threads $TASK_CPUS \
            --output-prefix "$P" \
            "${{INPUT[@]}}" \
            2>&1 | tee "$P.log" \
            || exit 1
        shopt -s nullglob
        mv ./*.tre ./*.rep ./*.one ./*.all ./*.unc ./*.log "$OUT"/
    """,
)

ganon_report_cmd = shell.create(
    name="ganon_report",
    image=GANON_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"rep": File, "db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        REP=({{inputs.rep}})
        DB={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {FIND_DB_PREFIX}
        W=$(mktemp -d); cd "$W"
        ganon \
            report \
            --input "${{REP[0]}}" \
            --output-prefix "$P" \
            --db-prefix "$DB_PREFIX" \
            $ARGS \
            || exit 1
        mv ./*.tre "$OUT"/
    """,
)

ganon_table_cmd = shell.create(
    name="ganon_table",
    image=GANON_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"tre": list[File], "tre_order": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        TRE=({inputs.tre})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra TRE_NAMES <<< {inputs.tre_order}; TRE_DIR=$(dirname "${TRE[0]}")
        TRE=(); for N in "${TRE_NAMES[@]}"; do TRE+=("$TRE_DIR/$N"); done
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        # ganon names each sample after its file, so link them under their names.
        NAMES=()
        for F in "${TRE[@]}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        ganon \
            table \
            --input "${NAMES[@]}" \
            --output-file "$P.txt" \
            $ARGS \
            || exit 1
        mv "$P.txt" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "ganon",
    ganon_buildcustom_cmd.as_task(),
    ganon_classify_cmd.as_task(),
    ganon_report_cmd.as_task(),
    ganon_table_cmd.as_task(),
)


@dataclass
class GanonBuildResult:
    filter: File  # <prefix>.hibf, or <prefix>.ibf
    tax: File | None  # <prefix>.tax
    info: File | None  # <prefix>.info.tsv, with --write-info-file
    db: Dir  # everything above


async def ganon_buildcustom(
    fasta: File | list[File],
    input_tsv: File | None = None,
    taxonomy_files: list[File] | None = None,
    genome_size_files: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> GanonBuildResult:
    """Build a ganon database of ``fasta``, listed by ``input_tsv`` if given."""
    fastas = fasta if isinstance(fasta, list) else [fasta]
    db = await ganon_buildcustom_cmd(
        fasta=fastas,
        input_tsv=input_tsv,
        taxonomy_files=taxonomy_files or [],
        genome_size_files=genome_size_files or [],
        prefix=prefix,
        args=args,
    )
    bloom = await db.get_file(f"{prefix}.hibf") or await db.get_file(f"{prefix}.ibf")
    if bloom is None:
        raise FileNotFoundError(f"ganon build-custom wrote no {prefix}.hibf or {prefix}.ibf")
    return GanonBuildResult(
        filter=bloom,
        tax=await db.get_file(f"{prefix}.tax"),
        info=await db.get_file(f"{prefix}.info.tsv"),
        db=db,
    )


@dataclass
class GanonClassifyResult:
    tre: File  # <prefix>.tre, the default taxonomic report
    rep: File  # <prefix>.rep
    log: File  # <prefix>.log
    one: File | None  # <prefix>.one, with --output-one
    all: File | None  # <prefix>.all, with --output-all
    unc: File | None  # <prefix>.unc, with --output-unclassified
    results: Dir  # everything above


async def ganon_classify(
    reads_1: File, db: Dir, reads_2: File | None = None, prefix: str = "sample", args: str = ""
) -> GanonClassifyResult:
    """Classify one sample's reads (``reads_2`` None for single-end) against the ganon database ``db``."""
    results = await ganon_classify_cmd(reads_1=reads_1, reads_2=reads_2, db=db, prefix=prefix, args=args)
    tre = await results.get_file(f"{prefix}.tre")
    rep = await results.get_file(f"{prefix}.rep")
    log = await results.get_file(f"{prefix}.log")
    if tre is None or rep is None or log is None:
        raise FileNotFoundError(f"ganon classify wrote no {prefix}.tre, {prefix}.rep or {prefix}.log")
    return GanonClassifyResult(
        tre=tre,
        rep=rep,
        log=log,
        one=await results.get_file(f"{prefix}.one"),
        all=await results.get_file(f"{prefix}.all"),
        unc=await results.get_file(f"{prefix}.unc"),
        results=results,
    )


async def ganon_report(rep: File, db: Dir, prefix: str = "sample", args: str = "") -> File:
    """A taxonomic report of ganon classify's ``rep`` against ``db``: ``<prefix>.tre``."""
    results = await ganon_report_cmd(rep=rep, db=db, prefix=prefix, args=args)
    tre = await results.get_file(f"{prefix}.tre")
    if tre is None:
        raise FileNotFoundError(f"ganon report wrote no {prefix}.tre")
    return tre


async def ganon_table(tre: File | list[File], prefix: str = "sample", args: str = "") -> File:
    """Several samples' ganon reports combined into one table: ``<prefix>.txt``."""
    reports = tre if isinstance(tre, list) else [tre]
    results = await ganon_table_cmd(tre=reports, tre_order=" ".join(f.name for f in reports), prefix=prefix, args=args)
    table = await results.get_file(f"{prefix}.txt")
    if table is None:
        raise FileNotFoundError(f"ganon table wrote no {prefix}.txt")
    return table
