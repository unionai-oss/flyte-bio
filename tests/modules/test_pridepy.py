"""Tests for flyte_bio.modules.pridepy.

The cases and expected md5s are upstream's pridepy/downloadfile and
pridepy/fetchsdrf tests (modules commit 9a1cc9f, pridepy 0.0.15). Both
download project PXD004684's SDRF from PRIDE when they run.
"""

from flyte.io import File

from flyte_bio.modules.pridepy import pridepy_downloadfile, pridepy_fetchsdrf
from tests.framework import assert_md5, env

SDRF_MD5 = "a8a2aa91d85adb241110fe3a7ce8f372"


async def check_sdrf(sdrf: File, name: str) -> None:
    assert sdrf.path.endswith("/" + name), sdrf.path
    async with sdrf.open("rb") as fh:
        assert "Source Name" in bytes(await fh.read()).decode(), "SDRF lacks 'Source Name'"
    await assert_md5(sdrf, SDRF_MD5, label=name)


@env.task
async def test_pridepy_downloadfile() -> None:
    # upstream case: pridepy/downloadfile "pride - PXD004684 - sdrf.tsv"
    await check_sdrf(await pridepy_downloadfile("sdrf.tsv", "PXD004684"), "sdrf.tsv")


@env.task
async def test_pridepy_fetchsdrf() -> None:
    # upstream case: pridepy/fetchsdrf "pride - PXD004684"
    await check_sdrf(await pridepy_fetchsdrf("PXD004684", prefix="PXD004684"), "PXD004684.sdrf.tsv")


tests = [test_pridepy_downloadfile, test_pridepy_fetchsdrf]
