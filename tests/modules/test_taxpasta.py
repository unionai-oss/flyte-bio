"""Tests for flyte_bio.modules.taxpasta.

The cases and expected md5s are upstream's taxpasta/standardise and
taxpasta/merge tests (modules commit 9a1cc9f, taxpasta 0.7.0): the test
data's Kraken2 report standardised, and merged with itself as two samples
(``test_1``, ``test_2``), without and with a samplesheet. As upstream, the
merged table's md5 is taken over its lines sorted, each line's fields
sorted, joined by newlines (Groovy's ``split`` drops trailing empty fields).
"""

import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.taxpasta import taxpasta_merge, taxpasta_standardise
from tests.framework import assert_md5, env, fixture

REPORT = "genomics/sarscov2/metagenome/test_1.kraken2.report.txt"


async def taxpasta_local_file(name: str, text: str) -> File:
    local = Path(tempfile.mkdtemp()) / name
    local.write_text(text)
    return await File.from_local(str(local))


async def taxpasta_test_reports() -> list[File]:
    async with (await fixture(REPORT)).open("rb") as fh:
        text = bytes(await fh.read()).decode()
    return [
        await taxpasta_local_file("test_1.kraken2.report.txt", text),
        await taxpasta_local_file("test_2.kraken2.report.txt", text),
    ]


def groovy_fields(line: str) -> list[str]:
    """``line.split("\\t")`` as Java splits it: trailing empty fields dropped."""
    fields = line.split("\t")
    while fields and fields[-1] == "":
        fields.pop()
    return fields


async def check_sorted_md5(table: File, expected: str) -> None:
    async with table.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    canonical = "\n".join("\t".join(sorted(groovy_fields(line))) for line in sorted(lines))
    actual = hashlib.md5(canonical.encode()).hexdigest()
    assert actual == expected, f"{table.name}: sorted md5 {actual}, expected {expected}"


@env.task
async def test_taxpasta_standardise() -> None:
    # upstream case: taxpasta/standardise "sarscov2 - metagenome - kraken report"
    table = await taxpasta_standardise(await fixture(REPORT), "kraken2", "tsv", prefix="test")
    assert table.name == "test.tsv", table.name
    await assert_md5(table, "6b62032ed403f431eeb0e17464ccc69b", label="test.tsv")


@env.task
async def test_taxpasta_merge() -> None:
    # upstream case: taxpasta/merge "sarscov2 - metagenome - kraken report"
    table = await taxpasta_merge(await taxpasta_test_reports(), "kraken2", "tsv", prefix="kraken2")
    assert table.name == "kraken2.tsv", table.name
    await check_sorted_md5(table, "e6b2f3988546fc2b56a69b1cb8fadfb9")


@env.task
async def test_taxpasta_merge_samplesheet() -> None:
    # upstream case: taxpasta/merge "sarscov2 - metagenome - kraken report - samplesheet"
    sheet = await taxpasta_local_file(
        "samplesheet.tsv",
        'sample\tprofile\ntest_1\t"test_1.kraken2.report.txt"\ntest_2\t"test_2.kraken2.report.txt"\n',
    )
    table = await taxpasta_merge(await taxpasta_test_reports(), "kraken2", "tsv", samplesheet=sheet, prefix="kraken2")
    await check_sorted_md5(table, "223fa8b2726c668a5df3e043067c019a")


tests = [test_taxpasta_standardise, test_taxpasta_merge, test_taxpasta_merge_samplesheet]
