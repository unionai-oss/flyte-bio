"""pycoqc — QC report of a nanopore run from its sequencing summary.

Exposes :func:`pycoqc`, which reads a basecaller's
``sequencing_summary.txt`` and writes an interactive HTML report,
``<prefix>.html``, and its statistics as JSON, ``<prefix>.json``. ``args``
carries pycoQC's options (e.g. ``--min_pass_qual 0``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PYCOQC_IMAGE = "quay.io/biocontainers/pycoqc:2.5.2--py_0"
# Upstream: process_medium; pycoQC is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

pycoqc_cmd = shell.create(
    name="pycoqc",
    image=PYCOQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"summary": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SUMMARY=({inputs.summary})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        pycoQC \
            $ARGS \
            -f "${SUMMARY[0]}" \
            -o "$P.html" \
            -j "$P.json" || exit 1
        mv "$P.html" "$P.json" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("pycoqc", pycoqc_cmd.as_task())


@dataclass
class PycoqcResult:
    html: File  # <prefix>.html
    json: File  # <prefix>.json


async def pycoqc(summary: File, prefix: str = "sample", args: str = "") -> PycoqcResult:
    """pycoQC's report on the sequencing summary ``summary``."""
    results = await pycoqc_cmd(summary=summary, prefix=prefix, args=args)
    picked = []
    for name in (f"{prefix}.html", f"{prefix}.json"):
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"pycoQC wrote no {name}")
        picked.append(found)
    return PycoqcResult(html=picked[0], json=picked[1])
