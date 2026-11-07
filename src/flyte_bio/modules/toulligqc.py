"""toulligqc — Oxford Nanopore run QC reports with ToulligQC.

Exposes :func:`toulligqc`, which reports on one sample's ONT data, passed to
ToulligQC by its extension as upstream: reads (``.fastq``, ``.fastq.gz``,
``.fq``, ``.fq.gz``) as ``--fastq``, a sequencing summary (``.txt``,
``.txt.gz``) as ``--sequencing-summary-source``, a BAM as ``--bam``. Any
other input is passed to none of them: give its option in ``args``. Pass
``--barcoding --barcodes ...`` in ``args`` for a barcoded run.

ToulligQC writes its report into a ``<prefix>/`` directory, kept as upstream:

- ``<prefix>/report.data``: the statistics, ``key=value`` per line;
- ``<prefix>/report.html``: the HTML report;
- ``<prefix>/images/*.html`` and ``<prefix>/images/plotly.min.js``: the
  interactive plots.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

TOULLIGQC_IMAGE = "community.wave.seqera.io/library/toulligqc_huggingface_hub_matplotlib_numpy_pruned:190ebbd12957a8dd"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

toulligqc_cmd = shell.create(
    name="toulligqc",
    image=TOULLIGQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"ont_file": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ONT=({inputs.ont_file})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        INPUT=()
        case "${ONT[0]}" in
            *.fastq|*.fastq.gz|*.fq|*.fq.gz) INPUT=(--fastq "${ONT[0]}") ;;
            *.txt|*.txt.gz) INPUT=(--sequencing-summary-source "${ONT[0]}") ;;
            *.bam) INPUT=(--bam "${ONT[0]}") ;;
        esac
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib"
        toulligqc \
            "${INPUT[@]}" \
            --output-directory . \
            --report-name "$P" \
            $ARGS
        mv "$P" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("toulligqc", toulligqc_cmd.as_task())


@dataclass
class ToulligqcResult:
    report_data: File  # <prefix>/report.data
    report_html: File | None  # <prefix>/report.html
    plots_html: list[File]  # <prefix>/images/*.html
    plotly_js: File  # <prefix>/images/plotly.min.js
    results: Dir  # the <prefix>/ directory and everything in it


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def toulligqc(ont_file: File, prefix: str = "sample", args: str = "") -> ToulligqcResult:
    """ToulligQC report on ``ont_file``: reads, a sequencing summary or a BAM."""
    results = await toulligqc_cmd(ont_file=ont_file, prefix=prefix, args=args)
    files = sorted([f async for f in results.walk()], key=lambda f: f.path)

    def in_report(f: File, sub: str = "") -> bool:
        return os.path.dirname(f.path).endswith(f"/{prefix}{sub}")

    data = [f for f in files if in_report(f) and file_name(f).endswith(".data")]
    html = [f for f in files if in_report(f) and file_name(f).endswith(".html")]
    plots = [f for f in files if in_report(f, "/images") and file_name(f).endswith(".html")]
    plotly = [f for f in files if in_report(f, "/images") and file_name(f) == "plotly.min.js"]
    if not data or not plotly:
        raise FileNotFoundError(f"toulligqc wrote no {prefix}/*.data or {prefix}/images/plotly.min.js")
    return ToulligqcResult(
        report_data=data[0],
        report_html=html[0] if html else None,
        plots_html=plots,
        plotly_js=plotly[0],
        results=results,
    )
