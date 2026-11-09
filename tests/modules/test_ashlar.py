"""Tests for flyte_bio.modules.ashlar.

The cases and expected md5s are upstream's ashlar tests (modules commit
9a1cc9f, ashlar 1.19.0): one and three CyCIF tonsil cycles, with and without
dark- and flat-field profiles. The md5s are over the output with its UUIDs
blanked, as the module does. As in upstream's snapshot, the profiles leave
these outputs unchanged.
"""

from flyte.io import File

from flyte_bio.modules.ashlar import ashlar
from tests.framework import assert_md5, env, fixture

OME_TIFF = "imaging/ome-tiff/"
ONE_CYCLE_MD5 = "bcf2234d791690e826583e71df3c7a75"
ALL_CYCLES_MD5 = "dd7617d7fa0da24c078cef185566ba99"


async def ashlar_cycles(n: int) -> list[File]:
    return [await fixture(f"{OME_TIFF}cycif-tonsil-cycle{i}.ome.tif") for i in range(1, n + 1)]


async def ashlar_profiles() -> tuple[list[File], list[File]]:
    return [await fixture(OME_TIFF + "cycif-tonsil-dfp.ome.tif")], [
        await fixture(OME_TIFF + "cycif-tonsil-ffp.ome.tif")
    ]


async def check(output: File, name: str, md5: str) -> None:
    assert output.name == name, f"got {output.name}, expected {name}"
    await assert_md5(output, md5, label=name)


@env.task
async def test_ashlar_with_fields() -> None:
    # upstream case: "ashlar - with fields"
    dfp, ffp = await ashlar_profiles()
    output = await ashlar(await ashlar_cycles(1), dfp=dfp, ffp=ffp, prefix="test_all")
    await check(output, "test_all.ome.tif", ONE_CYCLE_MD5)


@env.task
async def test_ashlar_one_cycle() -> None:
    # upstream case: "ashlar - 1 cycle - ome-tif"
    output = await ashlar(await ashlar_cycles(1), prefix="test_1")
    await check(output, "test_1.ome.tif", ONE_CYCLE_MD5)


@env.task
async def test_ashlar_all_cycles() -> None:
    # upstream case: "ashlar - all cycles - ome-tif"
    output = await ashlar(await ashlar_cycles(3), prefix="test_1")
    await check(output, "test_1.ome.tif", ALL_CYCLES_MD5)


@env.task
async def test_ashlar_all_cycles_with_fields() -> None:
    # upstream case: "ashlar - all cycles - with fields - ome-tif"
    dfp, ffp = await ashlar_profiles()
    output = await ashlar(await ashlar_cycles(3), dfp=dfp, ffp=ffp, prefix="test_1")
    await check(output, "test_1.ome.tif", ALL_CYCLES_MD5)


tests = [test_ashlar_with_fields, test_ashlar_one_cycle, test_ashlar_all_cycles, test_ashlar_all_cycles_with_fields]
