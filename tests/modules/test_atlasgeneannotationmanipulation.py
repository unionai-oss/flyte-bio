"""Tests for flyte_bio.modules.atlasgeneannotationmanipulation.

The cases and expected md5s are upstream's
atlasgeneannotationmanipulation/gtf2featureannotation tests (modules commit
9a1cc9f), with the FASTA case's test config. nf-test hashes the gzipped
cDNAs decompressed.
"""

from flyte_bio.modules.atlasgeneannotationmanipulation import atlasgeneannotationmanipulation_gtf2featureannotation
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"


@env.task
async def test_atlasgeneannotationmanipulation_gtf() -> None:
    # upstream case: "human - gtf"
    result = await atlasgeneannotationmanipulation_gtf2featureannotation(
        await fixture(GENOME + "genome.gtf"), prefix="homo_sapiens"
    )
    await assert_md5(result.feature_annotation, "9e4373ab72f445f948621a29d834f8a7", label="anno.tsv")
    assert result.filtered_cdna is None, "cDNAs written without a FASTA"


@env.task
async def test_atlasgeneannotationmanipulation_gtf_fasta() -> None:
    # upstream case: "human - gtf - fasta"
    result = await atlasgeneannotationmanipulation_gtf2featureannotation(
        await fixture(GENOME + "genome.gtf"),
        prefix="homo_sapiens_test",
        fasta=await fixture(GENOME + "transcriptome.fasta"),
        args="--feature-type transcript --filter-cdnas-output homo_sapiens_test.filtered.fa.gz",
    )
    await assert_md5(result.feature_annotation, "5159cc9ee2b3646b0e633b739aad97ff", label="anno.tsv")
    assert result.filtered_cdna is not None and result.filtered_cdna.name == "homo_sapiens_test.filtered.fa.gz"
    await assert_gunzipped_md5(result.filtered_cdna, "10c60e38188d7bbb697e703941652f04", label="filtered cDNAs")


tests = [test_atlasgeneannotationmanipulation_gtf, test_atlasgeneannotationmanipulation_gtf_fasta]
