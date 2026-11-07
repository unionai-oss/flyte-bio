"""sequali — sequencing-read QC reports.

Exposes :func:`sequali`, which reports on one sample's reads (a FASTQ, a
FASTQ pair, or an unaligned BAM) and writes ``<prefix>.html`` and
``<prefix>.json`` into an output Dir.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SEQUALI_IMAGE = "quay.io/biocontainers/sequali:1.0.2--py312h0fa9677_0"

# Upstream's process_low: two CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

sequali_cmd = shell.create(
    name="sequali",
    image=SEQUALI_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for r in "${{R1[@]}}" "${{R2[@]}}"; do ln -s "$r" .; NAMES+=("$(basename "$r")"); done
        sequali $ARGS -t $TASK_CPUS --html "$OUT/$P.html" --json "$OUT/$P.json" "${{NAMES[@]}}" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("sequali", sequali_cmd.as_task())


@dataclass
class SequaliResult:
    html: File  # <prefix>.html
    json: File  # <prefix>.json


async def sequali(reads_1: File, reads_2: File | None, prefix: str, args: str = "") -> SequaliResult:
    """Sequali's HTML and JSON QC reports for ``reads_1`` (and mate ``reads_2``)."""
    results = await sequali_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)
    html, json = await results.get_file(f"{prefix}.html"), await results.get_file(f"{prefix}.json")
    if html is None or json is None:
        raise FileNotFoundError(f"sequali wrote no {prefix}.html / {prefix}.json")
    return SequaliResult(html=html, json=json)
