"""Tests for flyte_bio.modules.knotannotsv.

The cases and expected values are upstream's knotannotsv tests (modules
commit 9a1cc9f, knotAnnotSV 1.1.5), with their params' ``--genomeBuild hg38``.
As upstream, they report the AnnotSV table fixture ``sv_query.tsv`` and check
the HTML report's line count and the workbook's size.
"""

from flyte.io import File

from flyte_bio.modules.knotannotsv import knotannotsv
from tests.framework import env, fixture

SV_TABLE = "genomics/homo_sapiens/illumina/tsv/sv_query.tsv"
ARGS = "--genomeBuild hg38"


async def knotannotsv_bytes(file: File) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read())


@env.task
async def test_knotannotsv_html() -> None:
    # upstream case: knotannotsv "homo_sapiens - annotsv_tsv, knot_out_xl - HTML output"
    report = await knotannotsv(await fixture(SV_TABLE), prefix="test", args=ARGS)
    assert report.path.endswith("/test.html"), report.path
    lines = len((await knotannotsv_bytes(report)).decode().splitlines())
    assert lines == 5249, f"test.html has {lines} lines"


@env.task
async def test_knotannotsv_xlsm() -> None:
    # upstream case: knotannotsv "homo_sapiens - annotsv_tsv, knot_out_xl - XLSM output"
    report = await knotannotsv(await fixture(SV_TABLE), prefix="test", excel=True, args=ARGS)
    assert report.path.endswith("/test.xlsm"), report.path
    size = len(await knotannotsv_bytes(report))
    assert size >= 89000, f"test.xlsm is {size} bytes"


tests = [test_knotannotsv_html, test_knotannotsv_xlsm]
