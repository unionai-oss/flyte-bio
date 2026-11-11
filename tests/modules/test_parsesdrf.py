"""Tests for flyte_bio.modules.parsesdrf.

The cases and expected md5s are upstream's parsesdrf/convert tests (modules
commit 9a1cc9f, sdrf-pipelines 0.1.4), all on the PXD009752 SDRF from the
mhcquant test data. As upstream, the MaxQuant parameter file isn't
reproducible (Python's hash seed reorders its modifications), so that case
checks the output names only.
"""

from flyte.io import File

from flyte_bio.modules.parsesdrf import parsesdrf_convert
from tests.framework import assert_md5, env, fixture

MHCQUANT_DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/mhcquant/testdata/"
SDRF = MHCQUANT_DATA + "PXD009752.sdrf.tsv"
FASTA = MHCQUANT_DATA + "UP000005640_9606_500prot.fasta"


async def check_outputs(files: list[File], expected: dict[str, str | None]) -> None:
    names = [f.path.rsplit("/", 1)[-1] for f in files]
    assert names == list(expected), f"outputs {names}"
    for file, (name, md5) in zip(files, expected.items(), strict=True):
        if md5 is not None:
            await assert_md5(file, md5, label=name)


@env.task
async def test_parsesdrf_convert_mhcquant() -> None:
    # upstream case: parsesdrf/convert "sdrf - mhcquant - PXD009752"
    files = await parsesdrf_convert(await fixture(SDRF), "mhcquant", prefix="PXD009752")
    await check_outputs(
        files,
        {
            "PXD009752_samplesheet.tsv": "4d7d2f1ca8d3def75b05ae8c01c481d8",
            "PXD009752_presets.tsv": "3ba3b5c24f71025bfc1799ea87cd0a12",
        },
    )


@env.task
async def test_parsesdrf_convert_openms() -> None:
    # upstream case: parsesdrf/convert "sdrf - openms - PXD009752"
    files = await parsesdrf_convert(await fixture(SDRF), "openms", prefix="PXD009752")
    await check_outputs(
        files,
        {
            "PXD009752_samplesheet.tsv": "9216c3ad7db89fedede8f65fdb27c3fc",
            "PXD009752_experimental_design.tsv": "86dc1f997ac2eb733ef4b3082063f5cf",
        },
    )


@env.task
async def test_parsesdrf_convert_maxquant() -> None:
    # upstream case: parsesdrf/convert "sdrf - maxquant - PXD009752"
    files = await parsesdrf_convert(
        await fixture(SDRF), "maxquant", fasta=await fixture(FASTA), raw_folder=".", prefix="PXD009752"
    )
    await check_outputs(files, {"PXD009752.xml": None, "PXD009752_design.txt": None})


@env.task
async def test_parsesdrf_convert_msstats() -> None:
    # upstream case: parsesdrf/convert "sdrf - msstats - PXD009752"
    files = await parsesdrf_convert(await fixture(SDRF), "msstats", prefix="PXD009752")
    await check_outputs(files, {"PXD009752.csv": "0619e7ce6d4b44534825e86b24534e4a"})


@env.task
async def test_parsesdrf_convert_normalyzerde() -> None:
    # upstream case: parsesdrf/convert "sdrf - normalyzerde - PXD009752"
    files = await parsesdrf_convert(await fixture(SDRF), "normalyzerde", prefix="PXD009752")
    await check_outputs(
        files,
        {
            "PXD009752_design.csv": "7a59adfe619191555051105b9c1feacc",
            "PXD009752_comparisons.csv": "81051bcc2cf1bedf378224b0a93e2877",
        },
    )


@env.task
async def test_parsesdrf_convert_diann() -> None:
    # upstream case: parsesdrf/convert "sdrf - diann - PXD009752"
    files = await parsesdrf_convert(await fixture(SDRF), "diann", prefix="PXD009752")
    await check_outputs(
        files,
        {
            "PXD009752.cfg": "bf4bc46a652861a436efeb22e0d1b746",
            "PXD009752_design.tsv": "6e4348aeaecaf5d05521f5d67d1df9aa",
        },
    )


tests = [
    test_parsesdrf_convert_mhcquant,
    test_parsesdrf_convert_openms,
    test_parsesdrf_convert_maxquant,
    test_parsesdrf_convert_msstats,
    test_parsesdrf_convert_normalyzerde,
    test_parsesdrf_convert_diann,
]
