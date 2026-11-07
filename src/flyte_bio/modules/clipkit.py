"""clipkit — trim multiple sequence alignments with ClipKIT.

Exposes :func:`clipkit`, which trims the alignment ``aln``. Outputs keep
upstream's names, where ``<ext>`` is ``out_format`` (default ``clipkit``; it
only names the files, ClipKIT's output format is ``-of`` in ``args``):

- the trimmed alignment: ``<prefix>.<ext>``;
- ClipKIT's per-site log (``--log``): ``<prefix>.<ext>.log``;
- ClipKIT's console output: ``<prefix>.<ext>.txt``;
- with ``--complementary`` in ``args``: the trimmed-away sites,
  ``<prefix>.<ext>.complement``;
- with ``--report_json`` in ``args``: ``<prefix>.<ext>.report.json``.

``auxiliary_file`` is ClipKIT's ``--auxiliary_file`` (the site list of
``--mode cst``). Unlike upstream, a ClipKIT failure fails the task (upstream
pipes it through ``tee`` without ``pipefail``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CLIPKIT_IMAGE = "quay.io/biocontainers/clipkit:2.11.4--pyhdfd78af_0"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

clipkit_cmd = shell.create(
    name="clipkit",
    image=CLIPKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"aln": File, "out_format": str, "auxiliary_file": File | None, "prefix": str, "args": str},
    defaults={"out_format": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        ALN=({{inputs.aln}})
        shopt -s nullglob; AUX=({{inputs.auxiliary_file}}); shopt -u nullglob
        P={{inputs.prefix}}
        EXT={{inputs.out_format}}
        EXT=${{EXT:-clipkit}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        AUXOPT=()
        if [ ${{#AUX[@]}} -gt 0 ]; then AUXOPT=(--auxiliary_file "${{AUX[0]}}"); fi
        clipkit \
            $ARGS \
            "${{AUXOPT[@]}}" \
            --threads $TASK_CPUS \
            "${{ALN[0]}}" \
            --output "$P.$EXT" \
            --log \
            2>&1 | tee "$P.$EXT.txt" || exit 1
        mv "$P.$EXT"* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("clipkit", clipkit_cmd.as_task())


@dataclass
class ClipkitResult:
    alignment: File  # <prefix>.<ext>, the trimmed alignment
    log: File  # <prefix>.<ext>.log, the per-site log
    metadata: File  # <prefix>.<ext>.txt, ClipKIT's console output
    complementary: File | None  # <prefix>.<ext>.complement, with --complementary
    json: File | None  # <prefix>.<ext>.report.json, with --report_json


async def clipkit(
    aln: File, out_format: str = "", auxiliary_file: File | None = None, prefix: str = "sample", args: str = ""
) -> ClipkitResult:
    """Trim the alignment ``aln``: ``<prefix>.<ext>`` with its log and console output."""
    results = await clipkit_cmd(aln=aln, out_format=out_format, auxiliary_file=auxiliary_file, prefix=prefix, args=args)
    name = f"{prefix}.{out_format or 'clipkit'}"

    async def pick(suffix: str) -> File:
        found = await results.get_file(name + suffix)
        if found is None:
            raise FileNotFoundError(f"clipkit wrote no {name}{suffix}")
        return found

    return ClipkitResult(
        alignment=await pick(""),
        log=await pick(".log"),
        metadata=await pick(".txt"),
        complementary=await results.get_file(f"{name}.complement"),
        json=await results.get_file(f"{name}.report.json"),
    )
