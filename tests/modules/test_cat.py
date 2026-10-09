"""Tests for flyte_bio.modules.cat.

Expected md5s are the upstream cat/fastq snapshot values, which nf-test
computes over the decompressed content. The cases join shards in a fixed
order, so a reordering bug changes the md5.
"""



from flyte_bio.modules.cat import cat_fastq
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_cat_fastq_single_end() -> None:
    # upstream case: cat/fastq "test_cat_fastq_single_end"
    shards = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    merged = await cat_fastq(shards)
    await assert_gunzipped_md5(merged, "ee314a9bd568d06617171b0c85f508da", label="single_end")


@env.task
async def test_cat_fastq_paired_end() -> None:
    # upstream case: cat/fastq "test_cat_fastq_paired_end" — R1 and R2 shards merged separately
    r1 = await cat_fastq([await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test2_1.fastq.gz")])
    r2 = await cat_fastq([await fixture(FASTQ + "test_2.fastq.gz"), await fixture(FASTQ + "test2_2.fastq.gz")])
    await assert_gunzipped_md5(r1, "3ad9406595fafec8172368f9cd0b6a22", label="paired_end R1")
    await assert_gunzipped_md5(r2, "a52cab0b840c7178b0ea83df1fdbe8d5", label="paired_end R2")


tests = [test_cat_fastq_single_end, test_cat_fastq_paired_end]
