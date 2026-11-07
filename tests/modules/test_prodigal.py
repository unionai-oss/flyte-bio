"""Tests for flyte_bio.modules.prodigal.

The cases and expected md5s are upstream's prodigal tests (modules commit
9a1cc9f, Prodigal 2.6.3). nf-test computes the md5s over the decompressed
outputs.
"""

from flyte_bio.modules.prodigal import ProdigalResult, prodigal
from tests.framework import assert_gunzipped_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"


async def check_common(result: ProdigalResult) -> None:
    await assert_gunzipped_md5(result.nucleotide_fasta, "1bc8a05bcb72a3c324f5e4ffaa716d3b", label="test.fna.gz")
    await assert_gunzipped_md5(result.amino_acid_fasta, "7168b854103f3586ccfdb71a44c389f7", label="test.faa.gz")
    await assert_gunzipped_md5(result.all_gene_annotations, "e6d6c50f0c39e5169f84ae3c90837fa9", label="test_all.txt.gz")


@env.task
async def test_prodigal_gff() -> None:
    result = await prodigal(await fixture(GENOME), output_format="gff", prefix="test")
    await assert_gunzipped_md5(result.gene_annotations, "612c2724c2891c63350f171f74165757", label="test.gff.gz")
    await check_common(result)


@env.task
async def test_prodigal_gbk() -> None:
    result = await prodigal(await fixture(GENOME), output_format="gbk", prefix="test")
    await assert_gunzipped_md5(result.gene_annotations, "188b3a0e3f78740ded7f3ec4d876cb4b", label="test.gbk.gz")
    await check_common(result)


tests = [test_prodigal_gff, test_prodigal_gbk]
