"""Tests for flyte_bio.modules.coreograph.

The case and expected values are upstream's coreograph test (modules commit
9a1cc9f, UNet Coreograph 2.2.9): as upstream, each output is checked by its
name and size.
"""

from flyte.io import File

from flyte_bio.modules.coreograph import coreograph
from tests.framework import env, fixture


async def coreograph_size(file: File) -> int:
    size = 0
    async with file.open("rb") as fh:
        while chunk := await fh.read(1 << 20):
            size += len(chunk)
    return size


async def check(files: list[File], expected: list[tuple[str, int]], label: str) -> None:
    got = [(f.name, await coreograph_size(f)) for f in files]
    assert got == expected, f"{label}: got {got}, expected {expected}"


@env.task
async def test_coreograph_tif() -> None:
    # upstream case: "coreograph - tif"
    result = await coreograph(await fixture("imaging/core_detection/single_core_dapi.tif"))
    await check(result.cores, [("1.tif", 13518174)], "cores")
    await check(result.masks, [("1_mask.tif", 67518)], "masks")
    await check([result.tma_map], [("TMA_MAP.tif", 23219)], "tma_map")
    await check([result.centroids], [("centroidsY-X.txt", 22)], "centroids")


tests = [test_coreograph_tif]
