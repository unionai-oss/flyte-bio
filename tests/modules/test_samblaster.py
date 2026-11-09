"""Tests for flyte_bio.modules.samblaster.

The case is upstream's samblaster test (modules commit 9a1cc9f, samblaster
0.1.26): the output BAM's read sequences' md5 (nft-bam's ``md5Reads``).
"""

from flyte_bio.modules.samblaster import samblaster
from tests.framework import assert_reads_md5, env, fixture


@env.task
async def test_samblaster_paired_end_umi_unsorted_bam() -> None:
    # upstream case: samblaster "homo_sapiens-test_paired_end_umi_unsorted_bam" (no ext.args: the test loads no config)
    bam = await fixture("genomics/homo_sapiens/illumina/bam/umi/test.paired_end.umi_unsorted.bam")
    out = await samblaster(bam, prefix="test")
    assert out.name == "test.bam", out.name
    await assert_reads_md5(out, "c1b74864f32583faf7d9bcd82217ff4c", label="samblaster")


tests = [test_samblaster_paired_end_umi_unsorted_bam]
