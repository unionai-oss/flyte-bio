"""Tests for flyte_bio.modules.raven.

The case and expected md5s are upstream's raven test (modules commit 9a1cc9f,
Raven 1.6.1): the Bacteroides fragilis Nanopore reads. Gzipped outputs are
compared decompressed, as nf-test hashes them. The module's 4 CPUs match
upstream's test CPUs.
"""

from flyte_bio.modules.raven import raven
from tests.framework import assert_gunzipped_md5, env, fixture


@env.task
async def test_raven_bacteroides_fragilis() -> None:
    # upstream case: "test-raven"
    reads = await fixture("genomics/prokaryotes/bacteroides_fragilis/nanopore/fastq/test.fastq.gz")
    result = await raven(reads, prefix="test")
    await assert_gunzipped_md5(result.fasta, "775cca993a9dbcc09f7be65d801a696b", label="fasta")
    await assert_gunzipped_md5(result.gfa, "5342a527f24e7b15c5bd38c42ffd0be6", label="gfa")


tests = [test_raven_bacteroides_fragilis]
