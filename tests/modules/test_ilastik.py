"""Tests for flyte_bio.modules.ilastik.

The cases and expected values are upstream's ilastik/pixelclassification and
ilastik/multicut tests (modules commit 9a1cc9f, ilastik 1.4.0). As upstream,
the probabilities export is checked by md5, the other exports by their size,
and the multicut segmentation by its existence.
"""

from flyte.io import File

from flyte_bio.modules.ilastik import ilastik_multicut, ilastik_pixelclassification
from tests.framework import assert_md5, assert_nonempty, env, fixture

IMAGE = "imaging/h5/plant_wga.h5"
PIXEL_PROJECT = "imaging/ilp/plant_wga.pixel_prob.ilp"
EXPORT_SIZE = 463922


async def file_size(file: File) -> int:
    size = 0
    async with file.open("rb") as fh:
        while chunk := await fh.read(1 << 20):
            size += len(chunk)
    return size


async def classify(export_source: str) -> File:
    output = await ilastik_pixelclassification(
        await fixture(IMAGE), await fixture(PIXEL_PROJECT), "h5", export_source, prefix="test"
    )
    assert output.name == "test.h5", f"got {output.name}"
    return output


async def check_size(export_source: str) -> None:
    size = await file_size(await classify(export_source))
    assert size == EXPORT_SIZE, f"{export_source}: test.h5 is {size} bytes, expected {EXPORT_SIZE}"


@env.task
async def test_ilastik_export_probabilities() -> None:
    await assert_md5(await classify("probabilities"), "a6331b68581ed037a2edcbc82fa862a6", label="test.h5")


@env.task
async def test_ilastik_export_simple_segmentation() -> None:
    await check_size("simple segmentation")


@env.task
async def test_ilastik_export_uncertainty() -> None:
    await check_size("uncertainty")


@env.task
async def test_ilastik_export_features() -> None:
    await check_size("features")


@env.task
async def test_ilastik_multicut_h5() -> None:
    mask = await ilastik_multicut(
        await fixture(IMAGE),
        await fixture("imaging/h5/plant_wga_probabilities.h5"),
        await fixture("imaging/ilp/plant_wga.multicut.ilp"),
        prefix="test",
    )
    assert mask.name == "test.tiff", f"got {mask.name}"
    await assert_nonempty(mask, label="test.tiff")


tests = [
    test_ilastik_export_probabilities,
    test_ilastik_export_simple_segmentation,
    test_ilastik_export_uncertainty,
    test_ilastik_export_features,
    test_ilastik_multicut_h5,
]
