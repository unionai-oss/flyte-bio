"""coreograph — find and crop the cores of a tissue microarray (TMA) image.

Exposes :func:`coreograph`, which detects the cores of a TMA image (typically
its DNA channel) with UNet Coreograph's U-Net and writes, as upstream:

- one cropped image per core, ``<n>.tif``;
- one mask per core, ``masks/<n>_mask.tif``;
- ``TMA_MAP.tif``, the TMA with the detected cores outlined and numbered;
- ``centroidsY-X.txt``, the cores' centroids (y, x), one per line.

The cores and masks come back in core-number order.
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

COREOGRAPH_IMAGE = "docker.io/labsyspharm/unetcoreograph:2.2.9"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

# TensorFlow's thread pools follow the task's CPUs.
coreograph_cmd = shell.create(
    name="coreograph",
    image=COREOGRAPH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGE=({{inputs.image}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W"
        export OMP_NUM_THREADS=$TASK_CPUS TF_NUM_INTRAOP_THREADS=$TASK_CPUS TF_NUM_INTEROP_THREADS=$TASK_CPUS
        python /app/UNetCoreograph.py --imagePath "${{IMAGE[0]}}" --outputPath . $ARGS || exit 1
        shopt -s nullglob
        CORES=(./*[0-9]*.tif)
        [ ${{#CORES[@]}} -gt 0 ] || {{ echo "coreograph found no cores" >&2; exit 1; }}
        mkdir "$OUT/cores"
        mv "${{CORES[@]}}" "$OUT/cores"/
        mv masks TMA_MAP.tif centroidsY-X.txt "$OUT"/ || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("coreograph", coreograph_cmd.as_task())


@dataclass
class CoreographResult:
    cores: list[File]  # <n>.tif, one per core
    masks: list[File]  # <n>_mask.tif, one per core
    tma_map: File  # TMA_MAP.tif
    centroids: File  # centroidsY-X.txt
    results: Dir  # cores/, masks/, TMA_MAP.tif and centroidsY-X.txt


def core_number(file: File) -> tuple[int, str]:
    found = re.search(r"\d+", file.name)
    return (int(found.group()) if found else 0, file.name)


async def coreograph(image: File, args: str = "") -> CoreographResult:
    """The cores of the TMA ``image``, their masks, the TMA map and the core centroids."""
    results = await coreograph_cmd(image=image, args=args)
    cores, masks = [], []
    async for file in results.walk():
        parent = file.path.rstrip("/").split("/")[-2]
        if parent == "cores":
            cores.append(file)
        elif parent == "masks":
            masks.append(file)
    tma_map = await results.get_file("TMA_MAP.tif")
    centroids = await results.get_file("centroidsY-X.txt")
    if tma_map is None or centroids is None:
        raise FileNotFoundError("coreograph wrote no TMA_MAP.tif or centroidsY-X.txt")
    return CoreographResult(
        cores=sorted(cores, key=core_number),
        masks=sorted(masks, key=core_number),
        tma_map=tma_map,
        centroids=centroids,
        results=results,
    )
