"""Tests for flyte_bio.modules.ucsc and bedtools genomecov (the bigWig steps).

Expected md5s are the upstream bedtools/genomecov, ucsc/bedclip and
ucsc/bedgraphtobigwig snapshot values.
"""

from flyte_bio.modules.bedtools import bedtools_genomecov
from flyte_bio.modules.ucsc import bedclip, bedgraphtobigwig
from tests.framework import assert_md5, env, fixture

SARSCOV2 = "genomics/sarscov2/"


@env.task
async def test_genomecov_bam_sorted() -> None:
    # upstream case: bedtools/genomecov "sarscov2 - no scale" (BAM, sorted, txt)
    bam = await fixture(SARSCOV2 + "illumina/bam/test.paired_end.bam")
    out = await bedtools_genomecov(bam, prefix="test.coverage", extension="txt", sort=True)
    await assert_md5(out, "66083198daca6c001d328ba9616e9b53", label="genomecov")


@env.task
async def test_bedclip() -> None:
    # upstream case: ucsc/bedclip "sarscov2"
    bedgraph = await fixture(SARSCOV2 + "illumina/bedgraph/test.bedgraph")
    sizes = await fixture(SARSCOV2 + "genome/genome.sizes")
    out = await bedclip(bedgraph, sizes, prefix="test.clip")
    await assert_md5(out, "e02395e1f7c593b3f79563067159ebc2", label="bedclip")


@env.task
async def test_bedgraphtobigwig() -> None:
    # upstream case: ucsc/bedgraphtobigwig "Should run without failures"
    bedgraph = await fixture(SARSCOV2 + "illumina/bedgraph/test.bedgraph")
    sizes = await fixture(SARSCOV2 + "genome/genome.sizes")
    out = await bedgraphtobigwig(bedgraph, sizes, prefix="test")
    await assert_md5(out, "910ecc7f57e3bbd5fac5a8edba4f615d", label="bigWig")


tests = [test_genomecov_bam_sorted, test_bedclip, test_bedgraphtobigwig]
