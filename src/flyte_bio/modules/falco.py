"""falco — read quality-control reports (a C++ drop-in for FastQC).

Exposes :func:`falco`, which writes falco's FastQC-style reports for one
sample's reads (one or more FASTQs) into an output Dir.

As upstream, the naming depends on how many files are given:

- one file (single-end or interleaved): ``<prefix>_fastqc_data.txt``,
  ``<prefix>_summary.txt`` and ``<prefix>_report.html``;
- several files: falco's own names, after each input file:
  ``<file name>_fastqc_data.txt``, ``<file name>_summary.txt`` and
  ``<file name>_fastqc_report.html``.

The reads keep their file names (they're symlinked into the working
directory), so the reports name them as upstream's do.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's main.nf names the quay.io biocontainer directly.
FALCO_IMAGE = "quay.io/biocontainers/falco:1.2.5--h077b44d_0"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

falco_cmd = shell.create(
    name="falco",
    image=FALCO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        READS=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        FILES=()
        for r in "${{READS[@]}}"; do ln -s "$r" "$(basename "$r")"; FILES+=("$(basename "$r")"); done
        if [ ${{#FILES[@]}} -eq 1 ]; then
            falco $ARGS --threads $TASK_CPUS "${{FILES[0]}}" \
                -D "${{P}}_fastqc_data.txt" -S "${{P}}_summary.txt" -R "${{P}}_report.html"
        else
            falco $ARGS --threads $TASK_CPUS "${{FILES[@]}}"
        fi
        mv ./*.txt ./*.html "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("falco", falco_cmd.as_task())


@dataclass
class FalcoResult:
    html: list[File]  # the HTML reports, sorted by name
    txt: list[File]  # the data and summary tables, sorted by name
    results: Dir  # everything falco wrote


async def falco(reads: list[File], prefix: str = "sample", args: str = "") -> FalcoResult:
    """falco's QC reports for ``reads`` (one FASTQ, or several, e.g. both mates)."""
    results = await falco_cmd(reads=reads, prefix=prefix, args=args)
    files = sorted(await results.list_files(), key=lambda f: os.path.basename(f.path))
    return FalcoResult(
        html=[f for f in files if f.path.endswith(".html")],
        txt=[f for f in files if f.path.endswith(".txt")],
        results=results,
    )
