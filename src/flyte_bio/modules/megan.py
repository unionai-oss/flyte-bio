"""megan — MEGAN (Community Edition) tools for MEGAN's RMA files.

Exposes :func:`megan_rma2info`, which reports on an RMA6 file, such as
MALT's (``rma2info``): the report asked for in ``args`` (e.g. ``-c2c
Taxonomy``, read counts per taxon) as ``<prefix>.txt.gz`` and, with
``megan_summary``, a MEGAN summary file ``<prefix>.megan`` (``-es``).

The image ships MEGAN Community Edition (GPL-3.0).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MEGAN_IMAGE = "quay.io/biocontainers/megan:6.25.9--h9ee0642_0"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

megan_rma2info_cmd = shell.create(
    name="megan_rma2info",
    image=MEGAN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"rma6": File, "megan_summary": bool, "prefix": str, "args": str},
    defaults={"megan_summary": False, "args": ""},
    outputs={"results": Dir},
    script=r"""
        RMA=({inputs.rma6})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        SUMMARY=()
        if [ {inputs.megan_summary} = true ]; then SUMMARY=(-es "$P.megan"); fi
        # MEGAN is Java and may exit 255, which the cluster takes for an interruption.
        rma2info \
            -i "${RMA[0]}" \
            -o "$P.txt.gz" \
            "${SUMMARY[@]}" \
            $ARGS \
            || exit 1
        shopt -s nullglob
        mv ./*.txt.gz ./*.megan "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("megan", megan_rma2info_cmd.as_task())


@dataclass
class MeganRma2infoResult:
    txt: File  # <prefix>.txt.gz
    megan_summary: File | None  # <prefix>.megan, with megan_summary
    results: Dir  # both of the above


async def megan_rma2info(
    rma6: File, megan_summary: bool = False, prefix: str = "sample", args: str = ""
) -> MeganRma2infoResult:
    """Report on the RMA6 file ``rma6`` (the report chosen in ``args``): ``<prefix>.txt.gz``."""
    results = await megan_rma2info_cmd(rma6=rma6, megan_summary=megan_summary, prefix=prefix, args=args)
    txt = await results.get_file(f"{prefix}.txt.gz")
    if txt is None:
        raise FileNotFoundError(f"rma2info wrote no {prefix}.txt.gz")
    return MeganRma2infoResult(txt=txt, megan_summary=await results.get_file(f"{prefix}.megan"), results=results)
