"""Tests for flyte_bio.modules.kraken2.

The cases and expected md5s are upstream's kraken2 tests (modules commit
9a1cc9f, kraken2 2.1.6). Upstream untars the test database in a setup step;
here :func:`fixture_dir` extracts it. Upstream snapshots only the report and
the per-read classifications; for the saved reads it checks their names, as
these tests do.
"""

from flyte.io import Dir, File

from flyte_bio.modules.kraken2 import kraken2
from tests.framework import assert_md5, assert_nonempty, env, fixture, fixture_dir

FASTQ = "genomics/sarscov2/illumina/fastq/"
REPORT_MD5 = "4227755fe40478b8d7dc8634b489761e"


async def database() -> Dir:
    return await fixture_dir("genomics/sarscov2/genome/db/kraken2.tar.gz", "kraken2")


async def check_names(files: list[File], names: list[str], label: str) -> None:
    assert [f.name for f in files] == names, f"{label}: got {[f.name for f in files]}, expected {names}"
    for file in files:
        await assert_nonempty(file, label=f"{label} {file.name}")


@env.task
async def test_kraken2_single_end() -> None:
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await kraken2(reads, await database(), prefix="test", save_output_fastqs=True)
    await assert_md5(result.report, REPORT_MD5, label="report")
    await check_names(result.classified_reads_fastq, ["test.classified.fastq.gz"], "classified")
    await check_names(result.unclassified_reads_fastq, ["test.unclassified.fastq.gz"], "unclassified")
    assert result.classified_reads_assignment is None


@env.task
async def test_kraken2_paired_end() -> None:
    reads_1 = await fixture(FASTQ + "test_1.fastq.gz")
    reads_2 = await fixture(FASTQ + "test_2.fastq.gz")
    result = await kraken2(reads_1, await database(), reads_2=reads_2, prefix="test", save_output_fastqs=True)
    await assert_md5(result.report, REPORT_MD5, label="report")
    await check_names(
        result.classified_reads_fastq, ["test.classified_1.fastq.gz", "test.classified_2.fastq.gz"], "classified"
    )
    await check_names(
        result.unclassified_reads_fastq,
        ["test.unclassified_1.fastq.gz", "test.unclassified_2.fastq.gz"],
        "unclassified",
    )


@env.task
async def test_kraken2_single_end_save_reads_assignment() -> None:
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await kraken2(reads, await database(), prefix="test", save_reads_assignment=True)
    await assert_md5(result.report, REPORT_MD5, label="report")
    assert result.classified_reads_assignment is not None
    await assert_md5(result.classified_reads_assignment, "e7a90531f0d8d777316515c36fe4cae0", label="classifiedreads")
    assert result.classified_reads_fastq == [] and result.unclassified_reads_fastq == []
    assert await result.results.get_file("test.classified.fastq.gz") is None, "wrote reads without save_output_fastqs"


tests = [
    test_kraken2_single_end,
    test_kraken2_paired_end,
    test_kraken2_single_end_save_reads_assignment,
]
