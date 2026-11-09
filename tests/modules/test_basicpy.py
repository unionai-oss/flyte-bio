"""Tests for flyte_bio.modules.basicpy.

The case and expected values are upstream's basicpy test (modules commit
9a1cc9f, BaSiCPy 1.2.0 in the 1.2.0-patch5 image). As upstream, the profiles
are checked by name and size: the dark-field profile's exact size, the
flat-field profile's rounded to the nearest 1000 bytes.
"""

from flyte.io import File

from flyte_bio.modules.basicpy import basicpy
from tests.framework import env, fixture


async def basicpy_size(file: File) -> int:
    size = 0
    async with file.open("rb") as fh:
        while chunk := await fh.read(1 << 20):
            size += len(chunk)
    return size


@env.task
async def test_basicpy_ome_tiff() -> None:
    # upstream case: "basicpy - OME-TIFF"
    result = await basicpy(await fixture("imaging/ome-tiff/cycif-tonsil-cycle1.ome.tif"), prefix="test")
    assert result.dfp.name == "test-dfp.ome.tif", f"dark-field profile is {result.dfp.name}"
    assert result.ffp.name == "test-ffp.ome.tif", f"flat-field profile is {result.ffp.name}"
    dfp_size = await basicpy_size(result.dfp)
    assert dfp_size == 2684, f"test-dfp.ome.tif is {dfp_size} bytes, expected 2684"
    ffp_size = await basicpy_size(result.ffp)
    rounded = int(ffp_size / 1000 + 0.5) * 1000  # as Groovy's Math.round
    assert rounded == 352000, f"test-ffp.ome.tif is {ffp_size} bytes, expected about 352000"


tests = [test_basicpy_ome_tiff]
