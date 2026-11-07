"""Tests for flyte_bio.modules.bracken.

The cases and expected md5s are upstream's bracken/bracken tests (modules
commit 9a1cc9f, bracken 3.1, which reports itself as 3.0.1). As upstream,
each case first classifies the reads with kraken2 against the test
Kraken2/Bracken database, then runs Bracken on that report with the test's
``ext.args``: ``-l S -t 10 -r 150``, or ``-l G -t 10 -r 150`` for the genus
case.
"""

from flyte.io import Dir, File

from flyte_bio.modules.bracken import bracken
from flyte_bio.modules.kraken2 import kraken2
from tests.framework import assert_md5, env, fixture, fixture_dir

FASTQ = "genomics/sarscov2/illumina/fastq/"
SPECIES_ARGS = "-l S -t 10 -r 150"
GENUS_ARGS = "-l G -t 10 -r 150"


async def database() -> Dir:
    return await fixture_dir("genomics/sarscov2/genome/db/kraken2_bracken.tar.gz", "kraken2_bracken")


async def kraken_report(db: Dir, paired: bool) -> File:
    reads_1 = await fixture(FASTQ + "test_1.fastq.gz")
    reads_2 = await fixture(FASTQ + "test_2.fastq.gz") if paired else None
    return (await kraken2(reads_1, db, reads_2=reads_2, prefix="test")).report


async def check(paired: bool, args: str, tsv_md5: str, report_md5: str) -> None:
    db = await database()
    result = await bracken(await kraken_report(db, paired), db, prefix="test", args=args)
    assert result.report.name == "test.tsv", result.report.name
    await assert_md5(result.report, tsv_md5, label="test.tsv")
    await assert_md5(result.kraken_style_report, report_md5, label="test.kraken2.report_bracken.txt")


@env.task
async def test_bracken_single_end() -> None:
    await check(False, SPECIES_ARGS, "4a21ae14ff8a0311d55f139af5247838", "ca0fbeedc4353b5fdd081688823a33df")


@env.task
async def test_bracken_paired_end_genus() -> None:
    await check(True, GENUS_ARGS, "f609b09d6edb5ebc1ea1435d1dd46cde", "2ce58814420a3690da1f08e10e8d3a30")


@env.task
async def test_bracken_paired_end() -> None:
    await check(True, SPECIES_ARGS, "4a21ae14ff8a0311d55f139af5247838", "ca0fbeedc4353b5fdd081688823a33df")


tests = [
    test_bracken_single_end,
    test_bracken_paired_end_genus,
    test_bracken_paired_end,
]
