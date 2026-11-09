"""Tests for flyte_bio.modules.cdhit.

The case and expected md5s are upstream's cdhit/cdhitest test (modules
commit 9a1cc9f, CD-HIT 4.8.1): the SARS-CoV-2 contigs clustered with the
default options.
"""

from flyte_bio.modules.cdhit import cdhit_cdhitest
from tests.framework import assert_md5, env, fixture


@env.task
async def test_cdhit_cdhitest_contigs() -> None:
    # upstream case: cdhit/cdhitest "test-cdhit-cdhitest"
    result = await cdhit_cdhitest(await fixture("genomics/sarscov2/illumina/fasta/contigs.fasta"), prefix="test")
    assert result.fasta.name == "test.fa", result.fasta.name
    assert result.clusters.name == "test.fa.clstr", result.clusters.name
    await assert_md5(result.fasta, "c96caedf8c949a2b871b83152d99c652", label="test.fa")
    await assert_md5(result.clusters, "941a0ee41ff76036ecbd6fa3ac01be7a", label="test.fa.clstr")


tests = [test_cdhit_cdhitest_contigs]
