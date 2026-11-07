"""deeptools — coverage matrices, heatmaps and fingerprint plots for ChIP/ATAC-seq.

Exposes:

- :func:`deeptools_computematrix` — score bigWigs over regions (BED), writing
  ``<prefix>.computeMatrix.mat.gz`` (the matrix plotHeatmap and plotProfile
  read) and ``<prefix>.computeMatrix.vals.mat.tab`` (the values as a table).
  ``args`` starts with the mode, e.g. ``scale-regions -b 1000``.
- :func:`deeptools_plotheatmap` — a heatmap of a computeMatrix matrix:
  ``<prefix>.plotHeatmap.pdf`` and ``<prefix>.plotHeatmap.mat.tab``.
- :func:`deeptools_plotfingerprint` — a fingerprint plot of BAMs' cumulative
  coverage: ``<prefix>.plotFingerprint.pdf``, the raw counts
  ``<prefix>.plotFingerprint.raw.txt`` and the quality metrics
  ``<prefix>.plotFingerprint.qcmetrics.txt``.

deepTools labels samples after their file names, so the inputs are linked into
the working directory under their own names, in the caller's order (the order
of the plot's samples). As upstream, ``--numberOfProcessors`` is the task's
CPUs, and computeMatrix records it in the matrix header. plotFingerprint adds
``--extendReads <fragment_size>`` for single-end reads when ``fragment_size``
is set (upstream's ``params.fragment_size``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DEEPTOOLS_IMAGE = "community.wave.seqera.io/library/deeptools:3.5.6--511a2b05870b9ec3"
# Upstream's tests run with 4 CPUs, which computeMatrix records in its matrix header.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
HEATMAP_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

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

env = flyte.TaskEnvironment.from_task(
    "deeptools",
    deeptools_computematrix_cmd.as_task(),
    deeptools_plotheatmap_cmd.as_task(),
    deeptools_plotfingerprint_cmd.as_task(),
)


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
