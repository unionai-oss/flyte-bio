"""Tests for flyte_bio.modules.vireo.

The case and expected values are upstream's vireo test (modules commit
9a1cc9f): cellsnp-lite's cell genotypes (``--genotype``) of the chr21
demultiplexing BAM, two donors, with ``--forceLearnGT`` (its test config).
Every output is compared by md5, gzipped ones decompressed as nf-test hashes
them.
"""

from flyte_bio.modules.cellsnp import cellsnp_modea
from flyte_bio.modules.samtools import samtools_index
from flyte_bio.modules.vireo import vireo
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

DATA = "genomics/homo_sapiens/demultiplexing/"


@env.task
async def test_vireo_demultiplexing() -> None:
    # upstream case: vireo "demultiplexing"
    bam = await fixture(DATA + "chr21.bam")
    donors = await fixture(DATA + "donor_genotype_chr21.vcf")
    cells = await cellsnp_modea(
        bam,
        await samtools_index(bam=bam),
        await fixture(DATA + "barcodes.tsv"),
        prefix="sample1",
        region_vcf=donors,
        args="--genotype",
    )
    assert cells.cells is not None, "cellsnp-lite wrote no cell VCF"
    result = await vireo(2, donors, prefix="test", cell_data=cells.cells, args="--forceLearnGT")
    assert result.genotype_vcf is not None and result.filtered_variants is not None, "no donor genotypes"
    await assert_md5(result.donor_ids, "f9036fcf2d801c59dbabba7fb2c8c14d", label="donor_ids.tsv")
    await assert_md5(result.filtered_variants, "840cbbeb9d1aed7c699f5f8634c6933f", label="filtered_variants.tsv")
    await assert_gunzipped_md5(result.genotype_vcf, "b5336b91f2ce207cd9bf4273e58b9776", label="GT_donors.vireo.vcf.gz")
    await assert_gunzipped_md5(result.prob_doublets, "3aff77cd983f24e7954c81430662c577", label="prob_doublet.tsv.gz")
    await assert_gunzipped_md5(result.prob_singlets, "ee76b1a2fc7ae4fdd8ff6764baf09d2b", label="prob_singlet.tsv.gz")
    await assert_md5(result.summary, "94b0ce40cc86bbe7fbacebf98e685ee7", label="summary.tsv")


tests = [test_vireo_demultiplexing]
