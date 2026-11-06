"""nanoplot — long-read QC plots and summary statistics (NanoPlot).

Exposes :func:`nanoplot`, which runs NanoPlot on one sample's long reads
(``.fastq.gz`` / ``.fq.gz``, passed as ``--fastq``) or on an ONT sequencing
summary (``.txt``, passed as ``--summary``). Any other input is passed to
neither: give its input option in ``args``, as upstream does.

The outputs keep NanoPlot's names: ``NanoPlot-report.html`` and the
per-plot ``*.html`` files, the optional ``*.png`` images, and the
``NanoStats.txt`` statistics (all prefixed when ``args`` sets ``-p``). As
upstream, NanoPlot's dated log isn't kept. NanoPlot's ``-t`` is the task's
CPUs.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

NANOPLOT_IMAGE = "quay.io/biocontainers/nanoplot:1.48.0--pyhdfd78af_1"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

nanoplot_cmd = shell.create(
    name="nanoplot",
    image=NANOPLOT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        INPUT=()
        case "${{READS[0]}}" in
            *.fastq.gz|*.fq.gz) INPUT=(--fastq "${{READS[0]}}") ;;
            *.txt) INPUT=(--summary "${{READS[0]}}") ;;
        esac
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib"
        NanoPlot $ARGS -t $TASK_CPUS "${{INPUT[@]}}"
        shopt -s nullglob
        mv ./*.html ./*.png ./*.txt "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("nanoplot", nanoplot_cmd.as_task())


@dataclass
class NanoplotResult:
    html: list[File]  # NanoPlot-report.html and the interactive plots
    png: list[File]  # static plot images, when NanoPlot writes them
    txt: list[File]  # NanoStats.txt
    results: Dir  # everything kept


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def nanoplot(reads: File, args: str = "") -> NanoplotResult:
    """NanoPlot QC of ``reads``: a ``.fastq.gz`` / ``.fq.gz`` or a sequencing summary ``.txt``."""
    results = await nanoplot_cmd(reads=reads, args=args)
    files = sorted(await results.list_files(), key=file_name)

    def with_suffix(suffix: str) -> list[File]:
        return [f for f in files if file_name(f).endswith(suffix)]

    return NanoplotResult(html=with_suffix(".html"), png=with_suffix(".png"), txt=with_suffix(".txt"), results=results)
