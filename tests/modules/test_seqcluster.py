"""Tests for flyte_bio.modules.seqcluster.

The case and expected md5 are upstream's seqcluster/collapse test (modules
commit 9a1cc9f, seqcluster 1.2.9), with upstream's ``ext.args`` and prefix.
nf-test hashes the ``.gz`` output decompressed.
"""

from flyte_bio.modules.seqcluster import seqcluster_collapse
from tests.framework import assert_gunzipped_md5, env, fixture

READS = "https://github.com/nf-core/test-datasets/raw/smrnaseq/testdata/trimmed/small_Clone1_N1.fastp.fastq.gz"


@env.task
async def test_seqcluster_collapse_human_fastq() -> None:
    fastq = await seqcluster_collapse(await fixture(READS), prefix="test_1_seqcluster", args="-m 1 --min_size 15")
    await assert_gunzipped_md5(fastq, "21c736de10e306f14ec296eaeb38ef45", label="test_1_seqcluster.fastq.gz")


tests = [test_seqcluster_collapse_human_fastq]
