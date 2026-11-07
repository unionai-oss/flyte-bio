"""Tests for flyte_bio.modules.khmer.

The case and expected md5s are upstream's khmer/uniquekmers test (modules
commit 9a1cc9f, khmer 3.0.0a3): 50-mers of the SARS-CoV-2 genome.
"""

from flyte_bio.modules.khmer import khmer_uniquekmers
from tests.framework import assert_md5, env, fixture


@env.task
async def test_khmer_uniquekmers() -> None:
    # upstream case: khmer/uniquekmers "sarscov2 - fasta"
    result = await khmer_uniquekmers(await fixture("genomics/sarscov2/genome/genome.fasta"), 50, prefix="test")
    await assert_md5(result.report, "ee489abd3b244dea3640649e1790d55e", label="test.report.txt")
    await assert_md5(result.kmers, "496ebf23653a01c7a42d743e47c19f65", label="test.kmers.txt")


tests = [test_khmer_uniquekmers]
