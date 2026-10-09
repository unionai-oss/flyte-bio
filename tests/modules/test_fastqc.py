"""Tests for flyte_bio.modules.fastqc.

FastQC reports embed the run date, so (as upstream) these check the report
names and that the HTML report describes the reads, not md5s.
"""

from flyte.io import Dir

from flyte_bio.modules.fastqc import fastqc
from tests.framework import env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
REPORT_MARKER = "<tr><td>File type</td><td>Conventional base calls</td></tr>"


async def assert_report(results: Dir, stem: str) -> None:
    html = await results.get_file(f"{stem}_fastqc.html")
    zipped = await results.get_file(f"{stem}_fastqc.zip")
    assert html is not None and zipped is not None, f"missing {stem}_fastqc.html/.zip"
    async with html.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    assert REPORT_MARKER in text, f"{stem}_fastqc.html doesn't look like a FastQC report"


@env.task
async def test_fastqc_single_end() -> None:
    # upstream case: fastqc "sarscov2 single-end [fastq]"
    results = await fastqc(reads_1=await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    await assert_report(results, "test")


@env.task
async def test_fastqc_paired_end() -> None:
    # upstream case: fastqc "sarscov2 paired-end [fastq]"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    results = await fastqc(reads_1=r1, reads_2=r2, prefix="test")
    await assert_report(results, "test_1")
    await assert_report(results, "test_2")


tests = [test_fastqc_single_end, test_fastqc_paired_end]
