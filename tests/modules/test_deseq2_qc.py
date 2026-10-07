"""Tests for flyte_bio.modules.deseq2_qc.

Mirrors upstream's local deseq2_qc test ("parse count data correctly"),
which only checks that the PCA and sample-distance tables are produced
(DESeq2 output isn't byte-stable); additionally checks the MultiQC tables
carry the relabelled headers.
"""

from flyte_bio.modules.deseq2_qc import deseq2_qc
from tests.framework import env, fixture

DATA = "https://ngi-igenomes.s3.amazonaws.com/testdata/nf-core/pipelines/rnaseq/3.15/deseq2qc/"


@env.task
async def test_deseq2_qc() -> None:
    # upstream case: deseq2_qc "parse count data correctly"
    counts = await fixture(DATA + "countFile.tsv")
    pca_header = await fixture(DATA + "deseq2_pca_header.txt")
    clustering_header = await fixture(DATA + "deseq2_clustering_header.txt")
    result = await deseq2_qc(counts, pca_header=pca_header, clustering_header=clustering_header)
    for name in ("deseq2.pca.vals.txt", "deseq2.sample.dists.txt"):
        assert await result.results.get_file(name) is not None, f"no {name}"
    async with result.pca_multiqc.open("rb") as fh:
        pca = bytes(await fh.read()).decode()
    assert "star_salmon_deseq2_pca" in pca and "STAR_SALMON DESeq2 PCA" in pca, pca[:400]
    async with result.dists_multiqc.open("rb") as fh:
        dists = bytes(await fh.read()).decode()
    assert "star_salmon_deseq2_clustering" in dists, dists[:400]


tests = [test_deseq2_qc]
