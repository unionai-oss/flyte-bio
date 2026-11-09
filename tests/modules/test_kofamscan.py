"""Tests for flyte_bio.modules.kofamscan.

The cases are upstream's kofamscan tests (modules commit 9a1cc9f, KofamScan
1.3.0): the SARS-CoV-2 proteome against the test-datasets' small KOfam
profile set and KO list, with the test config's prefix ``<id>.kofamscan``,
in the default (``detail``) and ``--format detail-tsv`` formats. Upstream
checks only that the table is written (no md5); these also check its header.
"""

from flyte.io import File

from flyte_bio.modules.kofamscan import kofamscan
from tests.framework import env, fixture

DB = "genomics/sarscov2/genome/db/kofamscan/"


async def check_table(table: File, name: str) -> None:
    assert table.name == name, f"output is {table.name}, expected {name}"
    async with table.open("rb") as fh:
        first = bytes(await fh.read()).decode().splitlines()[0]
    assert first.startswith("#") and "gene name" in first, f"{name} header is {first!r}"


async def run_kofamscan(args: str) -> File:
    return await kofamscan(
        await fixture("genomics/sarscov2/genome/proteome.fasta"),
        await fixture(DB + "profiles.tar.gz"),
        await fixture(DB + "ko_list.gz"),
        prefix="test.kofamscan",
        args=args,
    )


@env.task
async def test_kofamscan_txt() -> None:
    await check_table(await run_kofamscan(""), "test.kofamscan.txt")


@env.task
async def test_kofamscan_tsv() -> None:
    await check_table(await run_kofamscan("--format detail-tsv"), "test.kofamscan.tsv")


tests = [
    test_kofamscan_txt,
    test_kofamscan_tsv,
]
