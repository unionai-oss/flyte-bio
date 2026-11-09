"""backsub — channel-wise background subtraction of a registered multiplexed image.

Exposes :func:`backsub`, which subtracts, for each channel whose ``background``
column in the markers CSV names another marker, that marker's (exposure-scaled)
signal, and drops the channels marked ``remove``. It writes the corrected
pyramidal image ``<prefix>.ome.tif`` and the matching markers table
``<prefix>.csv``. backsub fails if a ``background`` marker isn't in the
``marker_name`` column.

As upstream, ``prefix`` must differ from the image's name.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BACKSUB_IMAGE = "community.wave.seqera.io/library/pip_backsub:1f25bcdea7241dd7"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

backsub_cmd = shell.create(
    name="backsub",
    image=BACKSUB_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "markers": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGE=({{inputs.image}})
        MARKERS=({{inputs.markers}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib" NUMBA_CACHE_DIR="$W/.numba"
        backsub \
            --input "${{IMAGE[0]}}" \
            --markers "${{MARKERS[0]}}" \
            --output "$P.ome.tif" \
            --marker-output "$P.csv" \
            $ARGS || exit 1
        mv "$P.ome.tif" "$P.csv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("backsub", backsub_cmd.as_task())


@dataclass
class BacksubResult:
    image: File  # <prefix>.ome.tif
    markers: File  # <prefix>.csv


async def backsub(image: File, markers: File, prefix: str = "sample", args: str = "") -> BacksubResult:
    """Background-subtract ``image`` per its ``markers`` CSV."""
    if image.name == f"{prefix}.ome.tif":
        raise ValueError("backsub: the input and output names are the same; choose another prefix")
    results = await backsub_cmd(image=image, markers=markers, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"backsub wrote no {name}")
        return found

    return BacksubResult(image=await pick(f"{prefix}.ome.tif"), markers=await pick(f"{prefix}.csv"))
