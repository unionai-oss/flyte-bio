"""Tests for flyte_bio.modules.racon.

The case and expected md5 are upstream's racon test (modules commit 9a1cc9f,
racon 1.4.20); nf-test hashes the decompressed consensus.
"""

from flyte_bio.modules.racon import racon
from tests.framework import assert_gunzipped_md5, env, fixture

FRAGILIS = "genomics/prokaryotes/bacteroides_fragilis/"


@env.task
async def test_racon_single_end() -> None:
    reads = await fixture(FRAGILIS + "nanopore/fastq/test.fastq.gz")
    assembly = await fixture(FRAGILIS + "genome/genome.fna.gz")
    paf = await fixture(FRAGILIS + "genome/genome.paf")
    consensus = await racon(reads, assembly, paf, prefix="test")
    assert consensus.path.endswith("test_assembly_consensus.fasta.gz"), consensus.path
    await assert_gunzipped_md5(consensus, "e610bf2e1df1701b88675e36a844baf9", label="consensus")


tests = [test_racon_single_end]
