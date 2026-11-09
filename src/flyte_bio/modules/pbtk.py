"""pbtk — PacBio BAM toolkit.

Exposes :func:`pbtk_pbmerge`, which merges PacBio BAMs (e.g. a sample's
chunks) into ``<prefix>.bam``, with its PacBio index ``<prefix>.bam.pbi``
unless ``args`` turns it off (``--no-pbi``). ``-j`` follows the task's CPUs;
``args`` go to ``pbmerge`` and aren't shell-quoted.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PBTK_IMAGE = "quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0"
# Upstream's label is process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# As upstream, the inputs are staged into input/ under their own names.
pbtk_pbmerge_cmd = shell.create(
    name="pbtk_pbmerge",
    image=PBTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bams": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAMS=({{inputs.bams}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); mkdir "$W/input"; cd "$W"
        NAMES=()
        for F in "${{BAMS[@]}}"; do NAMES+=("input/$(basename "$F")"); ln -s "$F" "${{NAMES[-1]}}"; done
        pbmerge -o "$P.bam" -j $TASK_CPUS $ARGS "${{NAMES[@]}}" || exit 1
        shopt -s nullglob
        mv "$P.bam" "$P".bam.pb[i] {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("pbtk", pbtk_pbmerge_cmd.as_task())


@dataclass
class PbmergeResult:
    bam: File  # <prefix>.bam
    pbi: File | None  # <prefix>.bam.pbi
    results: Dir  # everything above


async def pbtk_pbmerge(bams: list[File], prefix: str = "sample", args: str = "") -> PbmergeResult:
    """Merge the PacBio ``bams`` into ``<prefix>.bam`` (and its ``.pbi``)."""
    results = await pbtk_pbmerge_cmd(bams=bams, prefix=prefix, args=args)
    bam = await results.get_file(f"{prefix}.bam")
    if bam is None:
        raise FileNotFoundError(f"pbmerge wrote no {prefix}.bam")
    return PbmergeResult(bam=bam, pbi=await results.get_file(f"{prefix}.bam.pbi"), results=results)
