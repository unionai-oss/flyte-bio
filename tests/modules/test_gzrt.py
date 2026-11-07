"""Tests for flyte_bio.modules.gzrt.

The cases and expected md5s are upstream's gzrt tests (modules commit 9a1cc9f,
gzrt 0.9.1): a FASTQ truncated after 10 kb, alone and paired with an intact
mate. nf-test hashes the recovered ``.fastq.gz`` decompressed.
"""

from flyte_bio.modules.gzrt import gzrt
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/homo_sapiens/illumina/fastq/"
CORRUPTED = FASTQ + "test2_1_corrupted_10kb.fastq.gz"
RECOVERED_1_MD5 = "6bc7d0b6304039eadba7321d21e7b2c9"


@env.task
async def test_gzrt_single_end() -> None:
    # upstream case: gzrt "Run gzrt on test2_1_corrupted single-end - fastq.gz"
    result = await gzrt(await fixture(CORRUPTED), prefix="test-single_recovered")
    assert [f.name for f in result.reads] == ["test-single_recovered.fastq.gz"], [f.name for f in result.reads]
    await assert_gunzipped_md5(result.reads[0], RECOVERED_1_MD5, label="recovered reads")


@env.task
async def test_gzrt_paired_end() -> None:
    # upstream case: gzrt "Run gzrt on test2_1_corrupted and test2_2 paired-end - fastq.gz"
    result = await gzrt(
        await fixture(CORRUPTED), await fixture(FASTQ + "test2_2.fastq.gz"), prefix="test-paired_recovered"
    )
    names = [f.name for f in result.reads]
    assert names == ["test-paired_recovered_1.fastq.gz", "test-paired_recovered_2.fastq.gz"], names
    await assert_gunzipped_md5(result.reads[0], RECOVERED_1_MD5, label="recovered mate 1")
    await assert_gunzipped_md5(result.reads[1], "641919946a0e572ca7483f6e54476f3b", label="recovered mate 2")


tests = [test_gzrt_single_end, test_gzrt_paired_end]
