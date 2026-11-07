"""Tests for flyte_bio.modules.fastqutils.

The cases and expected values are upstream's fastqutils/info tests (modules
commit 9a1cc9f, fastq_utils 0.25.3). Upstream's two failure cases check that
the process fails and that its error report carries fastq_info's message;
here the task's failure is caught and its error checked for the same message.
"""

from flyte_bio.modules.fastqutils import fastqutils_info
from tests.framework import assert_md5, env, fixture

HS_FASTQ = "genomics/homo_sapiens/illumina/fastq/"


async def expect_failure(reads: str, message: str) -> None:
    fastq = await fixture(HS_FASTQ + reads)
    try:
        await fastqutils_info(fastq, prefix="test")
    except Exception as e:
        assert message in str(e), f"{reads}: the error lacks {message!r}: {e}"
        return
    raise AssertionError(f"fastq_info accepted the malformed {reads}")


@env.task
async def test_fastqutils_info_truncated() -> None:
    # upstream case: fastqutils/info "human - fastq - truncated"
    await expect_failure("test_truncated_clean.fastq", "file truncated")


@env.task
async def test_fastqutils_info_quality_mismatch() -> None:
    # upstream case: fastqutils/info "human - fastq - quality-mismatch"
    await expect_failure("test_quality_mismatch.fastq", "sequence and quality don't have the same length 52!=47")


@env.task
async def test_fastqutils_info_paired_end() -> None:
    # upstream case: fastqutils/info "sarscov2 - fastq.gz - pe"
    reads_1 = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    reads_2 = await fixture("genomics/sarscov2/illumina/fastq/test_2.fastq.gz")
    report = await fastqutils_info(reads_1, reads_2, prefix="test")
    assert report.name == "test.txt", report.name
    await assert_md5(report, "63adccb1bf3b3dc697010c44d6706060", label="report")


tests = [test_fastqutils_info_truncated, test_fastqutils_info_quality_mismatch, test_fastqutils_info_paired_end]
