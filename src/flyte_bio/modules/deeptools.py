"""deeptools — coverage matrices, heatmaps and fingerprint plots for ChIP/ATAC-seq.

Exposes:

- :func:`deeptools_computematrix` — score bigWigs over regions (BED), writing
  ``<prefix>.computeMatrix.mat.gz`` (the matrix plotHeatmap and plotProfile
  read) and ``<prefix>.computeMatrix.vals.mat.tab`` (the values as a table).
  ``args`` starts with the mode, e.g. ``scale-regions -b 1000``.
- :func:`deeptools_plotheatmap` — a heatmap of a computeMatrix matrix:
  ``<prefix>.plotHeatmap.pdf`` and ``<prefix>.plotHeatmap.mat.tab``.
- :func:`deeptools_plotprofile` — a profile plot of a computeMatrix matrix:
  ``<prefix>.plotProfile.pdf`` and its data ``<prefix>.plotProfile.tab``.
- :func:`deeptools_plotfingerprint` — a fingerprint plot of BAMs' cumulative
  coverage: ``<prefix>.plotFingerprint.pdf``, the raw counts
  ``<prefix>.plotFingerprint.raw.txt`` and the quality metrics
  ``<prefix>.plotFingerprint.qcmetrics.txt``.
- :func:`deeptools_multibamsummary` — read coverage of BAMs over genome-wide
  bins (``multiBamSummary bins``), labelled by ``labels`` if given:
  ``<prefix>.bamSummary.npz`` (prefix ``all_bam`` by default), the matrix
  plotCorrelation and plotPCA read.
- :func:`deeptools_plotcorrelation` — a correlation plot of a multiBamSummary
  matrix (``spearman`` and ``heatmap`` unless ``method``/``plot_type`` say
  otherwise): ``<prefix>.plotCorrelation.pdf`` and the correlation matrix
  ``<prefix>.plotCorrelation.mat.tab``.
- :func:`deeptools_plotpca` — a PCA plot of a multiBamSummary matrix:
  ``<prefix>.plotPCA.pdf`` and its data ``<prefix>.plotPCA.tab``.
- :func:`deeptools_bamcoverage` — a BAM's (or CRAM's) coverage track:
  ``<prefix>.bigWig``, or ``<prefix>.bedgraph`` when ``args`` has
  ``--outFileFormat bedgraph`` (or ``-of bedgraph``). As upstream, a CRAM is
  first converted to BAM against its reference, since deepTools can't read
  CRAM.

deepTools labels samples after their file names, so the inputs are linked into
the working directory under their own names, in the caller's order (the order
of the plot's samples). As upstream, ``--numberOfProcessors`` is the task's
CPUs, and computeMatrix records it in the matrix header. plotFingerprint adds
``--extendReads <fragment_size>`` for single-end reads when ``fragment_size``
is set (upstream's ``params.fragment_size``).

bamCoverage's upstream module pins an image with samtools, so it has its own
environment; :data:`env` aggregates both.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DEEPTOOLS_IMAGE = "community.wave.seqera.io/library/deeptools:3.5.6--511a2b05870b9ec3"
DEEPTOOLS_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/deeptools_samtools:22ef5f45259c4336"
# Upstream's tests run with 4 CPUs, which computeMatrix records in its matrix header.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
HEATMAP_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
# Upstream's process_low.
BAMCOVERAGE_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# Inputs are linked into a working directory under their own names; matplotlib
# needs a writable HOME and config directory.
WORKDIR_SETUP = r"""
        W=$(mktemp -d); cd "$W"
        mkdir -p .home/.matplotlib
        export HOME="$W/.home" MPLCONFIGDIR="$W/.home/.matplotlib"
"""

# A list[File] is staged under the files' names and globbed alphabetically, so
# each list's names ("/"-separated) carry the caller's order into the script.
LINK_IN_ORDER = r"""
        link_in_order() {
            local DIR=$1 NAMES=$2 NAME
            local -a LIST
            IFS=/ read -ra LIST <<< "$NAMES"
            for NAME in "${LIST[@]}"; do
                ln -s "$DIR/$NAME" "$NAME"
                LINKED+=("$NAME")
            done
        }
"""

deeptools_computematrix_cmd = shell.create(
    name="deeptools_computematrix",
    image=DEEPTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bigwigs": list[File],
        "bigwig_names": str,
        "regions": list[File],
        "region_names": str,
        "blacklist": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BW=({{inputs.bigwigs}})
        BED=({{inputs.regions}})
        shopt -s nullglob; BLACKLIST=({{inputs.blacklist}}); shopt -u nullglob
        BW_NAMES={{inputs.bigwig_names}}
        BED_NAMES={{inputs.region_names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WORKDIR_SETUP}
        {LINK_IN_ORDER}
        LINKED=(); link_in_order "$(dirname "${{BW[0]}}")" "$BW_NAMES"; SCORES=("${{LINKED[@]}}")
        LINKED=(); link_in_order "$(dirname "${{BED[0]}}")" "$BED_NAMES"; REGIONS=("${{LINKED[@]}}")
        EXTRA=()
        if [ ${{#BLACKLIST[@]}} -gt 0 ]; then
            NAME=$(basename "${{BLACKLIST[0]}}")
            ln -s "${{BLACKLIST[0]}}" "$NAME"
            EXTRA+=(--blackListFileName "$NAME")
        fi
        computeMatrix \
            $ARGS \
            --regionsFileName "${{REGIONS[@]}}" \
            --scoreFileName "${{SCORES[@]}}" \
            --outFileName "$P.computeMatrix.mat.gz" \
            --outFileNameMatrix "$P.computeMatrix.vals.mat.tab" \
            --numberOfProcessors $TASK_CPUS \
            "${{EXTRA[@]}}"
        mv ./*.mat.gz ./*.mat.tab "$OUT"/
    """,
)

deeptools_plotheatmap_cmd = shell.create(
    name="deeptools_plotheatmap",
    image=DEEPTOOLS_IMAGE,
    resources=HEATMAP_RESOURCES,
    inputs={"matrix": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        MAT=({{inputs.matrix}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WORKDIR_SETUP}
        NAME=$(basename "${{MAT[0]}}")
        ln -s "${{MAT[0]}}" "$NAME"
        plotHeatmap \
            $ARGS \
            --matrixFile "$NAME" \
            --outFileName "$P.plotHeatmap.pdf" \
            --outFileNameMatrix "$P.plotHeatmap.mat.tab"
        mv ./*.pdf ./*.tab "$OUT"/
    """,
)

deeptools_plotprofile_cmd = shell.create(
    name="deeptools_plotprofile",
    image=DEEPTOOLS_IMAGE,
    resources=HEATMAP_RESOURCES,
    inputs={"matrix": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        MAT=({{inputs.matrix}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WORKDIR_SETUP}
        NAME=$(basename "${{MAT[0]}}")
        ln -s "${{MAT[0]}}" "$NAME"
        plotProfile \
            $ARGS \
            --matrixFile "$NAME" \
            --outFileName "$P.plotProfile.pdf" \
            --outFileNameData "$P.plotProfile.tab" || exit 1
        mv ./*.pdf ./*.tab "$OUT"/
    """,
)

# Each BAM's index is linked beside it.
deeptools_plotfingerprint_cmd = shell.create(
    name="deeptools_plotfingerprint",
    image=DEEPTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bams": list[File],
        "bam_names": str,
        "bais": list[File],
        "bai_names": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bams}})
        BAI=({{inputs.bais}})
        BAM_NAMES={{inputs.bam_names}}
        BAI_NAMES={{inputs.bai_names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WORKDIR_SETUP}
        {LINK_IN_ORDER}
        LINKED=(); link_in_order "$(dirname "${{BAM[0]}}")" "$BAM_NAMES"; BAMS=("${{LINKED[@]}}")
        LINKED=(); link_in_order "$(dirname "${{BAI[0]}}")" "$BAI_NAMES"
        plotFingerprint \
            $ARGS \
            --bamfiles "${{BAMS[@]}}" \
            --plotFile "$P.plotFingerprint.pdf" \
            --outRawCounts "$P.plotFingerprint.raw.txt" \
            --outQualityMetrics "$P.plotFingerprint.qcmetrics.txt" \
            --numberOfProcessors $TASK_CPUS
        mv ./*.pdf ./*.raw.txt ./*.qcmetrics.txt "$OUT"/
    """,
)

# Each BAM's index is linked beside it. `labels` is space-separated.
deeptools_multibamsummary_cmd = shell.create(
    name="deeptools_multibamsummary",
    image=DEEPTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bams": list[File],
        "bam_names": str,
        "bais": list[File],
        "bai_names": str,
        "labels": str,
        "blacklist": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"labels": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bams}})
        BAI=({{inputs.bais}})
        shopt -s nullglob; BLACKLIST=({{inputs.blacklist}}); shopt -u nullglob
        BAM_NAMES={{inputs.bam_names}}
        BAI_NAMES={{inputs.bai_names}}
        LABELS={{inputs.labels}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WORKDIR_SETUP}
        {LINK_IN_ORDER}
        LINKED=(); link_in_order "$(dirname "${{BAM[0]}}")" "$BAM_NAMES"; BAMS=("${{LINKED[@]}}")
        LINKED=(); link_in_order "$(dirname "${{BAI[0]}}")" "$BAI_NAMES"
        EXTRA=()
        if [ -n "$LABELS" ]; then EXTRA+=(--labels $LABELS); fi
        BLACKLIST_ARGS=()
        if [ ${{#BLACKLIST[@]}} -gt 0 ]; then
            NAME=$(basename "${{BLACKLIST[0]}}")
            ln -s "${{BLACKLIST[0]}}" "$NAME"
            BLACKLIST_ARGS+=(--blackListFileName "$NAME")
        fi
        multiBamSummary bins \
            $ARGS \
            "${{EXTRA[@]}}" \
            --bamfiles "${{BAMS[@]}}" \
            --numberOfProcessors $TASK_CPUS \
            --outFileName "$P.bamSummary.npz" \
            "${{BLACKLIST_ARGS[@]}}" || exit 1
        mv "$P.bamSummary.npz" "$OUT"/
    """,
)

deeptools_plotcorrelation_cmd = shell.create(
    name="deeptools_plotcorrelation",
    image=DEEPTOOLS_IMAGE,
    resources=HEATMAP_RESOURCES,
    inputs={"matrix": File, "method": str, "plot_type": str, "prefix": str, "args": str},
    defaults={"method": "spearman", "plot_type": "heatmap", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        MAT=({{inputs.matrix}})
        METHOD={{inputs.method}}
        PLOT_TYPE={{inputs.plot_type}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WORKDIR_SETUP}
        NAME=$(basename "${{MAT[0]}}")
        ln -s "${{MAT[0]}}" "$NAME"
        plotCorrelation \
            $ARGS \
            --corData "$NAME" \
            --corMethod "$METHOD" \
            --whatToPlot "$PLOT_TYPE" \
            --plotFile "$P.plotCorrelation.pdf" \
            --outFileCorMatrix "$P.plotCorrelation.mat.tab" || exit 1
        mv ./*.pdf ./*.tab "$OUT"/
    """,
)

deeptools_plotpca_cmd = shell.create(
    name="deeptools_plotpca",
    image=DEEPTOOLS_IMAGE,
    resources=HEATMAP_RESOURCES,
    inputs={"matrix": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        MAT=({{inputs.matrix}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WORKDIR_SETUP}
        NAME=$(basename "${{MAT[0]}}")
        ln -s "${{MAT[0]}}" "$NAME"
        plotPCA \
            $ARGS \
            --corData "$NAME" \
            --plotFile "$P.plotPCA.pdf" \
            --outFileNameData "$P.plotPCA.tab" || exit 1
        mv ./*.pdf ./*.tab "$OUT"/
    """,
)

# The alignment and its index, and the reference and its .fai, are linked side
# by side under their own names. A CRAM is converted to `<name>.bam` first, as
# upstream (deepTools can't read CRAM).
deeptools_bamcoverage_cmd = shell.create(
    name="deeptools_bamcoverage",
    image=DEEPTOOLS_SAMTOOLS_IMAGE,
    resources=BAMCOVERAGE_RESOURCES,
    inputs={
        "alignment": File,
        "index": File | None,
        "fasta": File | None,
        "fai": File | None,
        "blacklist": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}})
        shopt -s nullglob
        IDX=({{inputs.index}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); BLACKLIST=({{inputs.blacklist}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        case " $ARGS " in
            *" --outFileFormat bedgraph "*|*" -of bedgraph "*) EXT=bedgraph ;;
            *) EXT=bigWig ;;
        esac
        {WORKDIR_SETUP}
        NAME=$(basename "${{IN[0]}}"); ln -s "${{IN[0]}}" "$NAME"
        if [ ${{#IDX[@]}} -gt 0 ]; then ln -s "${{IDX[0]}}" "$(basename "${{IDX[0]}}")"; fi
        if [ "${{NAME##*.}}" = cram ]; then
            FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
            FAI_ARGS=()
            if [ ${{#FAI[@]}} -gt 0 ]; then
                ln -s "${{FAI[0]}}" "$(basename "${{FAI[0]}}")"
                FAI_ARGS=(--fai-reference "$(basename "${{FAI[0]}}")")
            fi
            BAM="${{NAME%.cram}}.bam"
            samtools view -T "$FA_NAME" "$NAME" "${{FAI_ARGS[@]}}" -@ $TASK_CPUS -o "$BAM" || exit 1
            samtools index -b "$BAM" -@ $TASK_CPUS || exit 1
        else
            BAM=$NAME
        fi
        BLACKLIST_ARGS=()
        if [ ${{#BLACKLIST[@]}} -gt 0 ]; then
            BL_NAME=$(basename "${{BLACKLIST[0]}}"); ln -s "${{BLACKLIST[0]}}" "$BL_NAME"
            BLACKLIST_ARGS=(--blackListFileName "$BL_NAME")
        fi
        bamCoverage \
            --bam "$BAM" \
            $ARGS \
            --numberOfProcessors $TASK_CPUS \
            --outFileName "$P.$EXT" \
            "${{BLACKLIST_ARGS[@]}}" || exit 1
        mv "$P.$EXT" "$OUT"/
    """,
)

# One image per environment: bamCoverage's has samtools too, so `env` aggregates two.
deeptools_tools_env = flyte.TaskEnvironment.from_task(
    "deeptools_tools",
    deeptools_computematrix_cmd.as_task(),
    deeptools_plotheatmap_cmd.as_task(),
    deeptools_plotprofile_cmd.as_task(),
    deeptools_plotfingerprint_cmd.as_task(),
    deeptools_multibamsummary_cmd.as_task(),
    deeptools_plotcorrelation_cmd.as_task(),
    deeptools_plotpca_cmd.as_task(),
)
deeptools_bamcoverage_env = flyte.TaskEnvironment.from_task(
    "deeptools_bamcoverage", deeptools_bamcoverage_cmd.as_task()
)
env = flyte.TaskEnvironment(name="deeptools", depends_on=[deeptools_tools_env, deeptools_bamcoverage_env])


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


def ordered_names(files: list[File], label: str) -> str:
    """The files' names, "/"-separated in list order; they must be distinct and nonempty."""
    if not files:
        raise ValueError(f"deeptools needs at least one {label}")
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"deeptools {label} need distinct names, got {names}")
    return "/".join(names)


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"deeptools wrote no {name}")
    return found


@dataclass
class ComputeMatrixResult:
    matrix: File  # <prefix>.computeMatrix.mat.gz
    table: File  # <prefix>.computeMatrix.vals.mat.tab
    results: Dir


async def deeptools_computematrix(
    bigwigs: list[File],
    regions: list[File],
    blacklist: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> ComputeMatrixResult:
    """Score ``bigwigs`` (in order) over ``regions``, skipping ``blacklist`` regions if given.

    ``args`` gives computeMatrix's mode and options, e.g. ``scale-regions -b 1000``.
    """
    results = await deeptools_computematrix_cmd(
        bigwigs=bigwigs,
        bigwig_names=ordered_names(bigwigs, "bigWigs"),
        regions=regions,
        region_names=ordered_names(regions, "region files"),
        blacklist=blacklist,
        prefix=prefix,
        args=args,
    )
    return ComputeMatrixResult(
        matrix=await pick(results, f"{prefix}.computeMatrix.mat.gz"),
        table=await pick(results, f"{prefix}.computeMatrix.vals.mat.tab"),
        results=results,
    )


@dataclass
class PlotHeatmapResult:
    pdf: File  # <prefix>.plotHeatmap.pdf
    table: File  # <prefix>.plotHeatmap.mat.tab
    results: Dir


async def deeptools_plotheatmap(matrix: File, prefix: str = "sample", args: str = "") -> PlotHeatmapResult:
    """Plot a heatmap of a computeMatrix ``matrix`` (``.mat.gz``)."""
    results = await deeptools_plotheatmap_cmd(matrix=matrix, prefix=prefix, args=args)
    return PlotHeatmapResult(
        pdf=await pick(results, f"{prefix}.plotHeatmap.pdf"),
        table=await pick(results, f"{prefix}.plotHeatmap.mat.tab"),
        results=results,
    )


@dataclass
class PlotProfileResult:
    pdf: File  # <prefix>.plotProfile.pdf
    table: File  # <prefix>.plotProfile.tab
    results: Dir


async def deeptools_plotprofile(matrix: File, prefix: str = "sample", args: str = "") -> PlotProfileResult:
    """Plot a profile of a computeMatrix ``matrix`` (``.mat.gz``)."""
    results = await deeptools_plotprofile_cmd(matrix=matrix, prefix=prefix, args=args)
    return PlotProfileResult(
        pdf=await pick(results, f"{prefix}.plotProfile.pdf"),
        table=await pick(results, f"{prefix}.plotProfile.tab"),
        results=results,
    )


@dataclass
class PlotFingerprintResult:
    pdf: File  # <prefix>.plotFingerprint.pdf
    matrix: File  # <prefix>.plotFingerprint.raw.txt
    metrics: File  # <prefix>.plotFingerprint.qcmetrics.txt
    results: Dir


async def deeptools_plotfingerprint(
    bams: list[File],
    bais: list[File],
    prefix: str = "sample",
    single_end: bool = False,
    fragment_size: int = 0,
    args: str = "",
) -> PlotFingerprintResult:
    """Fingerprint plot of ``bams`` (in order), each indexed by its ``.bai`` in ``bais``.

    For single-end reads, a positive ``fragment_size`` extends reads to it.
    """
    if single_end and fragment_size > 0:
        args = f"{args} --extendReads {fragment_size}"
    results = await deeptools_plotfingerprint_cmd(
        bams=bams,
        bam_names=ordered_names(bams, "BAMs"),
        bais=bais,
        bai_names=ordered_names(bais, "BAM indexes"),
        prefix=prefix,
        args=args,
    )
    return PlotFingerprintResult(
        pdf=await pick(results, f"{prefix}.plotFingerprint.pdf"),
        matrix=await pick(results, f"{prefix}.plotFingerprint.raw.txt"),
        metrics=await pick(results, f"{prefix}.plotFingerprint.qcmetrics.txt"),
        results=results,
    )


async def deeptools_multibamsummary(
    bams: list[File],
    bais: list[File],
    labels: list[str] | None = None,
    blacklist: File | None = None,
    prefix: str = "all_bam",
    args: str = "",
) -> File:
    """Coverage of ``bams`` (in order, each indexed by its ``.bai`` in ``bais``) over genome-wide bins.

    ``labels`` name the samples (one per BAM, no whitespace); ``blacklist`` regions are skipped.
    Returns ``<prefix>.bamSummary.npz``.
    """
    labels = labels or []
    if labels and len(labels) != len(bams):
        raise ValueError(f"deeptools needs one label per BAM, got {len(labels)} for {len(bams)} BAMs")
    if any(not label or any(c.isspace() for c in label) for label in labels):
        raise ValueError(f"deeptools labels must be nonempty and have no whitespace, got {labels}")
    results = await deeptools_multibamsummary_cmd(
        bams=bams,
        bam_names=ordered_names(bams, "BAMs"),
        bais=bais,
        bai_names=ordered_names(bais, "BAM indexes"),
        labels=" ".join(labels),
        blacklist=blacklist,
        prefix=prefix,
        args=args,
    )
    return await pick(results, f"{prefix}.bamSummary.npz")


@dataclass
class PlotCorrelationResult:
    pdf: File  # <prefix>.plotCorrelation.pdf
    matrix: File  # <prefix>.plotCorrelation.mat.tab
    results: Dir


async def deeptools_plotcorrelation(
    matrix: File, method: str = "", plot_type: str = "", prefix: str = "sample", args: str = ""
) -> PlotCorrelationResult:
    """Correlation plot of a multiBamSummary ``matrix`` (``.npz``).

    ``method`` (``spearman`` if empty) and ``plot_type`` (``heatmap`` if empty) are plotCorrelation's
    ``--corMethod`` and ``--whatToPlot``.
    """
    results = await deeptools_plotcorrelation_cmd(
        matrix=matrix,
        method=method or "spearman",
        plot_type=plot_type or "heatmap",
        prefix=prefix,
        args=args,
    )
    return PlotCorrelationResult(
        pdf=await pick(results, f"{prefix}.plotCorrelation.pdf"),
        matrix=await pick(results, f"{prefix}.plotCorrelation.mat.tab"),
        results=results,
    )


@dataclass
class PlotPCAResult:
    pdf: File  # <prefix>.plotPCA.pdf
    table: File  # <prefix>.plotPCA.tab
    results: Dir


async def deeptools_plotpca(matrix: File, prefix: str = "sample", args: str = "") -> PlotPCAResult:
    """PCA plot of a multiBamSummary ``matrix`` (``.npz``)."""
    results = await deeptools_plotpca_cmd(matrix=matrix, prefix=prefix, args=args)
    return PlotPCAResult(
        pdf=await pick(results, f"{prefix}.plotPCA.pdf"),
        table=await pick(results, f"{prefix}.plotPCA.tab"),
        results=results,
    )


async def deeptools_bamcoverage(
    alignment: File,
    prefix: str,
    index: File | None = None,
    fasta: File | None = None,
    fai: File | None = None,
    blacklist: File | None = None,
    args: str = "",
) -> File:
    """The coverage track of ``alignment`` (a BAM, or a CRAM with ``fasta``), skipping ``blacklist`` regions.

    Returns ``<prefix>.bigWig``, or ``<prefix>.bedgraph`` when ``args`` has ``--outFileFormat bedgraph``
    (or ``-of bedgraph``).
    """
    if alignment.name.endswith(".cram") and fasta is None:
        raise ValueError("deeptools bamCoverage needs the reference fasta for a CRAM")
    padded = f" {args} "
    extension = "bedgraph" if " --outFileFormat bedgraph " in padded or " -of bedgraph " in padded else "bigWig"
    results = await deeptools_bamcoverage_cmd(
        alignment=alignment,
        index=index,
        fasta=fasta,
        fai=fai,
        blacklist=blacklist,
        prefix=prefix,
        args=args,
    )
    return await pick(results, f"{prefix}.{extension}")
