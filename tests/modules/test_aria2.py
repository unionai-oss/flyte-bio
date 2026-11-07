"""Tests for flyte_bio.modules.aria2.

The case and expected md5 are upstream's aria2 test (modules commit 9a1cc9f,
aria2 1.37.0): download a BAM from the test datasets.
"""

from flyte_bio.modules.aria2 import aria2
from tests.framework import TEST_DATA_BASE, assert_md5, env


@env.task
async def test_aria2_sarscov2_single_end_bam() -> None:
    bam = await aria2(TEST_DATA_BASE + "genomics/sarscov2/illumina/bam/test.single_end.bam")
    assert bam.name == "test.single_end.bam", f"downloaded {bam.name}"
    await assert_md5(bam, "21afed4c3e007de5e007cc5cbaebede7")


tests = [
    test_aria2_sarscov2_single_end_bam,
]
