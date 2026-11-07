"""Tests for flyte_bio.modules.ucsc and bedtools genomecov (the bigWig steps).

Expected md5s are the upstream bedtools/genomecov, ucsc/bedclip,
ucsc/bedgraphtobigwig, ucsc/wigtobigwig and ucsc/gtftogenepred snapshot
values; the gtftogenepred case uses its test config's ``-genePredExt
-geneNameAsName2`` (modules commit 9a1cc9f, UCSC 482).
"""

from flyte_bio.modules.bedtools import bedtools_genomecov
from flyte_bio.modules.ucsc import bedclip, bedgraphtobigwig, gtftogenepred, wigtobigwig
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


@env.task
async def test_ucsc_wigtobigwig() -> None:
    # upstream case: ucsc/wigtobigwig "test-ucsc-wigtobigwig"
    wig = await fixture(SARSCOV2 + "illumina/wig/test.wig.gz")
    sizes = await fixture(SARSCOV2 + "genome/genome.sizes")
    out = await wigtobigwig(wig, sizes, prefix="test")
    assert out.name == "test.bw", out.name
    await assert_md5(out, "b64af7003665dc51fae958216b06ed95", label="wigToBigWig")


@env.task
async def test_ucsc_gtftogenepred_gtf() -> None:
    # upstream case: ucsc/gtftogenepred "test-ucsc-gtftogenepred"
    gtf = await fixture(SARSCOV2 + "genome/genome.gtf")
    result = await gtftogenepred(gtf, prefix="test", args="-genePredExt -geneNameAsName2")
    assert result.genepred.name == "test.genepred", result.genepred.name
    assert result.refflat is not None and result.refflat.name == "test.refflat", result.refflat
    await assert_md5(result.genepred, "779e4749efaf38da3443ddfde30cc76c", label="genePred")
    await assert_md5(result.refflat, "4101802f41d4cf7ee2667587da11bf42", label="refFlat")


tests = [
    test_genomecov_bam_sorted,
    test_bedclip,
    test_bedgraphtobigwig,
    test_ucsc_wigtobigwig,
    test_ucsc_gtftogenepred_gtf,
]
