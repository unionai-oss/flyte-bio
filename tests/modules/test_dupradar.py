"""Tests for flyte_bio.modules.dupradar.

Expected md5s are the upstream dupradar snapshot values for the duplication
matrix and the two MultiQC tables (the PDFs aren't byte-stable, so — as
upstream — only their presence is checked).
"""

from flyte_bio.modules.dupradar import DupradarResult, dupradar
from tests.framework import assert_md5, env, fixture

SARSCOV2 = "genomics/sarscov2/"
# As upstream's test config: the sarscov2 GTF has CDS but no exon features.
FEATURE_TYPE = "--feature_type CDS"


async def check(result: DupradarResult, matrix: str, intercept: str, curve: str, label: str) -> None:
    await assert_md5(result.dup_matrix, matrix, label=f"{label} dupMatrix")
    await assert_md5(result.intercept_mqc, intercept, label=f"{label} dup_intercept_mqc")
    await assert_md5(result.curve_mqc, curve, label=f"{label} duprateExpDensCurve_mqc")
    for name in ("test_duprateExpBoxplot.pdf", "test_expressionHist.pdf", "test_duprateExpDens.pdf"):
        assert await result.results.get_file(name) is not None, f"{label}: no {name}"


@env.task
async def test_dupradar_single_end() -> None:
    # upstream case: dupradar "sarscov2 - bam - single_end"
    bam = await fixture(SARSCOV2 + "illumina/bam/test.single_end.bam")
    gtf = await fixture(SARSCOV2 + "genome/genome.gtf")
    result = await dupradar(bam, gtf, prefix="test", strandedness="forward", single_end=True, args=FEATURE_TYPE)
    await check(
        result,
        "2beda4c140548a2b8c91bf6bde01ddc6",
        "602a30c09b79533abbbb76efbb168e3e",
        "b41f82f7d515d5e2b3e7993add90d26b",
        "single_end",
    )


@env.task
async def test_dupradar_paired_end() -> None:
    # upstream case: dupradar "sarscov2 - bam - paired_end"
    bam = await fixture(SARSCOV2 + "illumina/bam/test.paired_end.bam")
    gtf = await fixture(SARSCOV2 + "genome/genome.gtf")
    result = await dupradar(bam, gtf, prefix="test", strandedness="forward", single_end=False, args=FEATURE_TYPE)
    await check(
        result,
        "5a327feaba56e5ea96a2cb7e8577e196",
        "4071a87d8db9731e9a74e51f00523056",
        "af58f18dd239447f04fa4a313124a549",
        "paired_end",
    )


tests = [test_dupradar_single_end, test_dupradar_paired_end]
