"""Tests for flyte_bio.modules.cellsnp.

The case and expected values are upstream's cellsnp/modea test (modules
commit 9a1cc9f): the chr21 demultiplexing BAM, indexed, at the donors' SNPs.
As upstream, the barcodes and depth matrices are compared by md5 and the base
VCF checked for presence.
"""

from flyte_bio.modules.cellsnp import cellsnp_modea
from flyte_bio.modules.samtools import samtools_index
from tests.framework import assert_md5, assert_nonempty, env, fixture

DATA = "genomics/homo_sapiens/demultiplexing/"


@env.task
async def test_cellsnp_modea_genotyping() -> None:
    # upstream case: cellsnp/modea "genotyping"
    bam = await fixture(DATA + "chr21.bam")
    result = await cellsnp_modea(
        bam,
        await samtools_index(bam=bam),
        await fixture(DATA + "barcodes.tsv"),
        prefix="sample1",
        region_vcf=await fixture(DATA + "donor_genotype_chr21.vcf"),
    )
    await assert_nonempty(result.base, label="base.vcf.gz")
    await assert_md5(result.samples, "9e488782c1bcd63c37ee3d1c4c0a9217", label="samples.tsv")
    await assert_md5(result.allele_depth, "dbe2d13dca2717749554d1f2a8b85650", label="tag.AD.mtx")
    await assert_md5(result.depth_coverage, "29fbb6241c6c7b0a0fa31021e622c415", label="tag.DP.mtx")
    await assert_md5(result.depth_other, "1e3429950c59edec58a80a9b4ecda552", label="tag.OTH.mtx")


tests = [test_cellsnp_modea_genotyping]
