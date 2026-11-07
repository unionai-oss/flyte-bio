"""Tests for flyte_bio.modules.deeptools.

The cases and expected values are upstream's deeptools computematrix,
plotfingerprint and plotheatmap tests (modules commit 9a1cc9f, deepTools
3.5.6). nf-test hashes the ``.mat.gz`` matrix decompressed; its header records
the 4 processors upstream's tests run with, which the module's 4 CPUs match.
Upstream checks only the plots' names (PDFs aren't reproducible), and for
plotFingerprint whether a line of the metrics is exactly ``0.24184576629880325``
(it isn't).
"""

from flyte.io import File

from flyte_bio.modules.deeptools import (
    deeptools_computematrix,
    deeptools_plotfingerprint,
    deeptools_plotheatmap,
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


tests = [
    test_deeptools_computematrix_no_blacklist,
    test_deeptools_computematrix_blacklist,
    test_deeptools_plotfingerprint_bam,
    test_deeptools_plotheatmap_mat,
]
