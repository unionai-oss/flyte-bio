"""Tests for flyte_bio.modules.deeptools.

The cases and expected values are upstream's deeptools computematrix,
plotfingerprint and plotheatmap tests (modules commit 9a1cc9f, deepTools
3.5.6). nf-test hashes the ``.mat.gz`` matrix decompressed; its header records
the 4 processors upstream's tests run with, which the module's 4 CPUs match.
Upstream checks only the plots' names (PDFs aren't reproducible), and for
plotFingerprint whether a line of the metrics is exactly ``0.24184576629880325``
(it isn't).

The multibamsummary, plotcorrelation and plotpca cases are upstream's at the
same commit: two labelled human BAMs summarised over bins, then plotted.
Upstream checks the npz's and plots' names, and the plotCorrelation matrix's
md5 (the same for the default and explicit method and plot type).
"""

from flyte.io import File

from flyte_bio.modules.deeptools import (
    deeptools_computematrix,
    deeptools_multibamsummary,
    deeptools_plotcorrelation,
    deeptools_plotfingerprint,
    deeptools_plotheatmap,
    deeptools_plotpca,
)
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, env, fixture

BIGWIG = "genomics/sarscov2/illumina/bigwig/test.bigwig"
BED = "genomics/sarscov2/genome/bed/test.bed"
BLACKLIST = "genomics/homo_sapiens/genome/genome.blacklist_intervals.bed"
COMPUTEMATRIX_ARGS = "scale-regions -b 1000"  # upstream's tests/nextflow.config
MATRIX_MD5 = "2bc77c2c0bdc61463c3aaebdb48e4bad"
TABLE_MD5 = "19e22051cc44edb7db3e0f8345330d90"
BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def head(file: File, size: int) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read(size))


@env.task
async def test_deeptools_computematrix_no_blacklist() -> None:
    result = await deeptools_computematrix(
        [await fixture(BIGWIG)], [await fixture(BED)], prefix="test", args=COMPUTEMATRIX_ARGS
    )
    assert result.matrix.name == "test.computeMatrix.mat.gz", result.matrix.name
    await assert_gunzipped_md5(result.matrix, MATRIX_MD5, label="matrix")
    await assert_md5(result.table, TABLE_MD5, label="table")


@env.task
async def test_deeptools_computematrix_blacklist() -> None:
    result = await deeptools_computematrix(
        [await fixture(BIGWIG)],
        [await fixture(BED)],
        blacklist=await fixture(BLACKLIST),
        prefix="test",
        args=COMPUTEMATRIX_ARGS,
    )
    await assert_gunzipped_md5(result.matrix, MATRIX_MD5, label="matrix")
    await assert_md5(result.table, TABLE_MD5, label="table")


@env.task
async def test_deeptools_plotfingerprint_bam() -> None:
    result = await deeptools_plotfingerprint(
        [await fixture(BAM)], [await fixture(BAM + ".bai")], prefix="test", single_end=False
    )
    assert result.pdf.name == "test.plotFingerprint.pdf", result.pdf.name
    assert result.matrix.name == "test.plotFingerprint.raw.txt", result.matrix.name
    assert await head(result.pdf, 4) == b"%PDF", "plotFingerprint pdf isn't a PDF"
    raw = (await text(result.matrix)).splitlines()
    assert raw[:2] == ["#plotFingerprint --outRawCounts", "'test.paired_end.sorted.bam'"], raw[:2]
    metrics = (await text(result.metrics)).splitlines()
    assert "0.24184576629880325" not in metrics, "upstream's snapshot has no metrics line 0.24184576629880325"
    assert metrics[0].startswith("Sample\tAUC\t"), f"unexpected metrics header {metrics[0]!r}"
    assert len(metrics) == 2 and metrics[1].startswith("test.paired_end.sorted.bam\t"), metrics


@env.task
async def test_deeptools_plotheatmap_mat() -> None:
    matrix = await fixture("genomics/sarscov2/illumina/deeptools/test.computeMatrix.mat.gz")
    result = await deeptools_plotheatmap(matrix, prefix="test")
    assert result.pdf.name == "test.plotHeatmap.pdf", result.pdf.name
    assert result.table.name == "test.plotHeatmap.mat.tab", result.table.name
    assert await head(result.pdf, 4) == b"%PDF", "plotHeatmap pdf isn't a PDF"
    await assert_nonempty(result.table, label="plotHeatmap table")


BAM2 = "genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam"
CORRELATION_MD5 = "24e46b74bd580e22c48687a1843eaa5f"


async def bam_summary() -> File:
    """Upstream's multiBamSummary run: two BAMs labelled test_bam1 and test_bam2, no blacklist."""
    bams = [await fixture(BAM), await fixture(BAM2)]
    bais = [await fixture(BAM + ".bai"), await fixture(BAM2 + ".bai")]
    return await deeptools_multibamsummary(bams, bais, labels=["test_bam1", "test_bam2"])


@env.task
async def test_deeptools_multibamsummary_bam() -> None:
    # upstream case: deeptools/multibamsummary "homo_sampiens - bam"
    matrix = await bam_summary()
    assert matrix.name == "all_bam.bamSummary.npz", matrix.name
    assert await head(matrix, 2) == b"PK", "bamSummary isn't an npz (zip) archive"


async def check_plotcorrelation(method: str, plot_type: str) -> None:
    result = await deeptools_plotcorrelation(await bam_summary(), method=method, plot_type=plot_type, prefix="test")
    assert result.pdf.name == "test.plotCorrelation.pdf", result.pdf.name
    assert result.matrix.name == "test.plotCorrelation.mat.tab", result.matrix.name
    await assert_md5(result.matrix, CORRELATION_MD5, label="correlation matrix")


@env.task
async def test_deeptools_plotcorrelation_bam() -> None:
    # upstream case: deeptools/plotcorrelation "homo_sampiens - bam"
    await check_plotcorrelation("spearman", "heatmap")


@env.task
async def test_deeptools_plotcorrelation_no_method() -> None:
    # upstream case: deeptools/plotcorrelation "homo_sampiens - bam - no_method"
    await check_plotcorrelation("", "heatmap")


@env.task
async def test_deeptools_plotcorrelation_no_plot_type() -> None:
    # upstream case: deeptools/plotcorrelation "homo_sampiens - bam - no_plot_type"
    await check_plotcorrelation("spearman", "")


@env.task
async def test_deeptools_plotpca_bam() -> None:
    # upstream case: deeptools/plotpca "homo_sampiens - bam"
    result = await deeptools_plotpca(await bam_summary(), prefix="test")
    assert result.pdf.name == "test.plotPCA.pdf", result.pdf.name
    assert result.table.name == "test.plotPCA.tab", result.table.name
    assert await head(result.pdf, 4) == b"%PDF", "plotPCA pdf isn't a PDF"
    await assert_nonempty(result.table, label="plotPCA table")


tests = [
    test_deeptools_computematrix_no_blacklist,
    test_deeptools_computematrix_blacklist,
    test_deeptools_plotfingerprint_bam,
    test_deeptools_plotheatmap_mat,
    test_deeptools_multibamsummary_bam,
    test_deeptools_plotcorrelation_bam,
    test_deeptools_plotcorrelation_no_method,
    test_deeptools_plotcorrelation_no_plot_type,
    test_deeptools_plotpca_bam,
]
