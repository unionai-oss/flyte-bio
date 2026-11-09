"""quilt — QUILT imputation of low-coverage sequencing against a reference panel.

Exposes :func:`quilt_quilt`, which imputes (and phases) the samples in BAM or
CRAM files over one region of a chromosome, from a reference panel's IMPUTE
haplotypes and legend files, into ``<prefix>.vcf.gz`` and its ``.tbi``.
Optional inputs are a BAM list (``--bamlist``; otherwise the reads files are
listed in the given order), sample names (``--sampleNames_file``), a positions
file (``--posfile``), truth haplotypes (``--phasefile``) or genotypes
(``--genfile``) for accuracy reports, and a genetic map
(``--genetic_map_file``). CRAMs need the reference FASTA (``--reference``).

Outputs keep upstream's names: ``<prefix>.vcf.gz``, ``<prefix>.vcf.gz.tbi``,
and the ``RData/`` (e.g. with ``--save_prepared_reference=TRUE``) and
``plots/`` (with ``--make_plots=TRUE``) directories. QUILT's log (its stderr)
is kept as ``<prefix>.quilt.log``, which upstream leaves in the task's
``.command.err``. As upstream, ``--seed=1`` is added unless ``args`` sets a
seed; ``--nCores`` follows the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

QUILT_IMAGE = "community.wave.seqera.io/library/r-quilt_r-base:32b4dc10c1885bdf"
# Upstream labels QUILT process_single; more than one core can change its output.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

quilt_quilt_cmd = shell.create(
    name="quilt_quilt",
    image=QUILT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": list[File],
        "indexes": list[File],
        "list_option": str,
        "bamlist": File | None,
        "sample_names": File | None,
        "haplotypes": File,
        "legend": File,
        "posfile": File | None,
        "phasefile": File | None,
        "genfile": File | None,
        "chromosome": str,
        "region_start": str,
        "region_end": str,
        "ngen": str,
        "buffer": str,
        "genetic_map": File | None,
        "fasta": File | None,
        "fasta_index": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"region_start": "", "region_end": "", "buffer": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        HAP=({{inputs.haplotypes}}); LEGEND=({{inputs.legend}})
        shopt -s nullglob
        READS=({{inputs.reads}}); IDX=({{inputs.indexes}}); BAMLIST=({{inputs.bamlist}})
        NAMES_FILE=({{inputs.sample_names}}); POS=({{inputs.posfile}}); PHASE=({{inputs.phasefile}})
        GEN=({{inputs.genfile}}); M=({{inputs.genetic_map}}); FA=({{inputs.fasta}}); FAI=({{inputs.fasta_index}})
        shopt -u nullglob
        LIST_OPTION={{inputs.list_option}}
        CHR={{inputs.chromosome}}; START={{inputs.region_start}}; END={{inputs.region_end}}
        NGEN={{inputs.ngen}}; BUFFER={{inputs.buffer}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        link() {{ ln -s "$1" "$(basename "$1")"; basename "$1"; }}
        READ_NAMES=()
        for F in "${{READS[@]}}"; do READ_NAMES+=("$(link "$F")"); done
        for F in "${{IDX[@]}}"; do link "$F" > /dev/null; done
        if [ ${{#BAMLIST[@]}} -gt 0 ]; then
            BAM_LIST=$(link "${{BAMLIST[0]}}")
        else
            printf "%s\n" "${{READ_NAMES[@]}}" > all_files.txt
            BAM_LIST=all_files.txt
        fi
        EXTRA=()
        if [ "$LIST_OPTION" = cram ]; then
            FA_NAME=$(link "${{FA[0]}}")
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
            EXTRA+=(--reference="$FA_NAME" --cramlist="$BAM_LIST")
        else
            EXTRA+=(--bamlist="$BAM_LIST")
        fi
        if [ ${{#M[@]}} -gt 0 ]; then EXTRA+=(--genetic_map_file="$(link "${{M[0]}}")"); fi
        if [ ${{#POS[@]}} -gt 0 ]; then EXTRA+=(--posfile="$(link "${{POS[0]}}")"); fi
        if [ ${{#PHASE[@]}} -gt 0 ]; then EXTRA+=(--phasefile="$(link "${{PHASE[0]}}")"); fi
        if [ ${{#GEN[@]}} -gt 0 ]; then EXTRA+=(--genfile="$(link "${{GEN[0]}}")"); fi
        if [ ${{#NAMES_FILE[@]}} -gt 0 ]; then EXTRA+=(--sampleNames_file="$(link "${{NAMES_FILE[0]}}")"); fi
        EXTRA+=(--chr="$CHR")
        if [ -n "$START" ]; then EXTRA+=(--regionStart="$START"); fi
        if [ -n "$END" ]; then EXTRA+=(--regionEnd="$END"); fi
        if [ -n "$BUFFER" ]; then EXTRA+=(--buffer="$BUFFER"); fi
        if [[ "$ARGS" != *--seed* ]]; then ARGS="$ARGS --seed=1"; fi

        # Fontconfig needs a writable cache directory.
        export XDG_CACHE_HOME=$PWD/.cache
        QUILT.R \
            "${{EXTRA[@]}}" \
            --nGen="$NGEN" \
            --nCores=$TASK_CPUS \
            --outputdir="." \
            --reference_haplotype_file="$(link "${{HAP[0]}}")" \
            --reference_legend_file="$(link "${{LEGEND[0]}}")" \
            --output_filename="$P.vcf.gz" \
            $ARGS 2> "$P.quilt.log" || {{ tail -n 40 "$P.quilt.log" >&2; exit 1; }}
        mv "$P.vcf.gz" "$P.quilt.log" "$OUT"/
        for F in "$P.vcf.gz.tbi" RData plots; do if [ -e "$F" ]; then mv "$F" "$OUT"/; fi; done
    """,
)

env = flyte.TaskEnvironment.from_task("quilt", quilt_quilt_cmd.as_task())


@dataclass
class QuiltResult:
    vcf: File  # <prefix>.vcf.gz: the imputed genotypes
    tbi: File | None  # <prefix>.vcf.gz.tbi
    log: File  # <prefix>.quilt.log: QUILT's stderr
    rdata: list[File]  # RData/*: e.g. the prepared reference, with --save_prepared_reference=TRUE
    plots: list[File]  # plots/*, with --make_plots=TRUE
    results: Dir


async def quilt_quilt(
    reads: list[File],
    indexes: list[File],
    haplotypes: File,
    legend: File,
    chromosome: str,
    ngen: int,
    prefix: str,
    region_start: int | None = None,
    region_end: int | None = None,
    buffer: int | None = None,
    bamlist: File | None = None,
    sample_names: File | None = None,
    posfile: File | None = None,
    phasefile: File | None = None,
    genfile: File | None = None,
    genetic_map: File | None = None,
    fasta: File | None = None,
    fasta_index: File | None = None,
    args: str = "",
) -> QuiltResult:
    """Impute the samples in ``reads`` (BAMs or CRAMs, with their ``indexes``) over a region of ``chromosome``.

    ``haplotypes`` and ``legend`` are the reference panel's IMPUTE files;
    ``ngen`` is the number of generations since admixture (``--nGen``).
    CRAMs need ``fasta`` (and its ``fasta_index``).
    """
    extensions = sorted({f.name.rsplit(".", 1)[-1] for f in reads})
    if extensions not in (["bam"], ["cram"]):
        raise ValueError(f"QUILT needs all BAM or all CRAM reads, got {extensions}")
    if extensions == ["cram"] and fasta is None:
        raise ValueError("QUILT needs the reference FASTA to read CRAMs")
    results = await quilt_quilt_cmd(
        reads=reads,
        indexes=indexes,
        list_option=extensions[0],
        bamlist=bamlist,
        sample_names=sample_names,
        haplotypes=haplotypes,
        legend=legend,
        posfile=posfile,
        phasefile=phasefile,
        genfile=genfile,
        chromosome=chromosome,
        region_start="" if region_start is None else str(region_start),
        region_end="" if region_end is None else str(region_end),
        ngen=str(ngen),
        buffer="" if buffer is None else str(buffer),
        genetic_map=genetic_map,
        fasta=fasta,
        fasta_index=fasta_index,
        prefix=prefix,
        args=args,
    )
    vcf, log = await results.get_file(f"{prefix}.vcf.gz"), await results.get_file(f"{prefix}.quilt.log")
    if vcf is None or log is None:
        raise FileNotFoundError(f"QUILT wrote no {prefix}.vcf.gz")
    files = [f async for f in results.walk()]

    def under(directory: str) -> list[File]:
        return sorted((f for f in files if f"/{directory}/" in f.path), key=lambda f: f.name)

    return QuiltResult(
        vcf=vcf,
        tbi=await results.get_file(f"{prefix}.vcf.gz.tbi"),
        log=log,
        rdata=under("RData"),
        plots=under("plots"),
        results=results,
    )
