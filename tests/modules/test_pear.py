"""Tests for flyte_bio.modules.pear.

The case and expected md5s are upstream's pear test (modules commit 9a1cc9f,
PEAR 0.9.6); nf-test hashes the decompressed reads. Upstream checks only the
discarded reads' name, since the file is empty.
"""

from flyte_bio.modules.pear import pear
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_pear_paired_end() -> None:
    result = await pear(
        await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz"), prefix="test"
    )
    await assert_gunzipped_md5(result.assembled, "c66f0d68268d800a0ca9cd285fe62453", label="assembled")
    forward, reverse = result.unassembled
    await assert_gunzipped_md5(forward, "066b58ee61df35bd2b063a511ffd5d85", label="unassembled forward")
    await assert_gunzipped_md5(reverse, "ed35907306bdbd96d59144ddb16c951a", label="unassembled reverse")
    assert result.discarded.path.endswith("test.discarded.fastq.gz"), result.discarded.path


tests = [test_pear_paired_end]
