"""Tests for flyte_bio.modules.souporcell.

The case and expected value are upstream's souporcell test (modules commit
9a1cc9f): the chr21 demultiplexing BAM and barcodes, one cluster, against
chr21. As upstream, ``clusters.tsv`` is compared by md5 and the genotypes VCF
and ambient RNA estimate (not stable) checked for presence.
"""

from flyte_bio.modules.souporcell import souporcell
from tests.framework import assert_md5, assert_nonempty, env, fixture

DATA = "genomics/homo_sapiens/demultiplexing/"


@env.task
async def test_souporcell_bam_barcodes() -> None:
    # upstream case: souporcell "homo sapiens - bam and barcodes"
    result = await souporcell(
        await fixture(DATA + "chr21.bam"),
        await fixture(DATA + "barcodes.tsv"),
        await fixture("genomics/homo_sapiens/genome/chr21/sequence/genome.fasta"),
        clusters=1,
        prefix="test",
    )
    await assert_md5(result.clusters, "dd4199a3dcb1cf949b96ab3497189426", label="clusters.tsv")
    await assert_nonempty(result.vcf, label="cluster_genotypes.vcf")
    await assert_nonempty(result.ambient_rna, label="ambient_rna.txt")


tests = [test_souporcell_bam_barcodes]
