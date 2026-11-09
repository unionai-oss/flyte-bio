"""basicpy — estimate illumination (shading) profiles of a tiled microscopy image.

Exposes :func:`basicpy`, which fits BaSiC's dark-field and flat-field
profiles to the tiles of one image (e.g. one imaging cycle), writing
``<prefix>-dfp.ome.tif`` and ``<prefix>-ffp.ome.tif``. They feed
:func:`flyte_bio.modules.ashlar.ashlar`'s illumination correction.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BASICPY_IMAGE = "docker.io/labsyspharm/basicpy-docker-mcmicro:1.2.0-patch5"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

basicpy_cmd = shell.create(
    name="basicpy",
    image=BASICPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGE=({{inputs.image}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib" NUMBA_CACHE_DIR="$W/.numba" XDG_CACHE_HOME="$W/.cache"
        /opt/main.py -i "${{IMAGE[0]}}" -o . --output-flatfield "$P" --output-darkfield "$P" $ARGS || exit 1
        mv "$P-dfp.ome.tif" "$P-ffp.ome.tif" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("basicpy", basicpy_cmd.as_task())


@dataclass
class BasicpyResult:
    dfp: File  # <prefix>-dfp.ome.tif, the dark-field profile
    ffp: File  # <prefix>-ffp.ome.tif, the flat-field profile


async def basicpy(image: File, prefix: str = "sample", args: str = "") -> BasicpyResult:
    """Dark- and flat-field profiles of ``image``'s tiles."""
    results = await basicpy_cmd(image=image, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"basicpy wrote no {name}")
        return found

    return BasicpyResult(dfp=await pick(f"{prefix}-dfp.ome.tif"), ffp=await pick(f"{prefix}-ffp.ome.tif"))
