"""Tests for flyte_bio.modules.mifaser.

The case is upstream's mifaser test (modules commit 9a1cc9f, mi-faser 1.64):
single-end reads against the ``GS-24-all`` database on upstream's
funcprofiler test-data branch, passed as its archive, with the three tables
checked by md5.
"""

from flyte_bio.modules.mifaser import mifaser
from tests.framework import assert_md5, env, fixture

DATABASE = "https://github.com/nf-core/test-datasets/raw/refs/heads/funcprofiler/data/database/mifaser/GS-24-all.tar.gz"


@env.task
async def test_mifaser_single_end_fastq() -> None:
    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    result = await mifaser([reads], await fixture(DATABASE), prefix="test")
    await assert_md5(result.multi_ec, "08ad0cac1771a190e63cd76a3b3ab686", label="multi_ec")
    await assert_md5(result.analysis, "722399e3bfd7af4273b36779a111db8a", label="analysis")
    await assert_md5(result.ec_counts, "195e3e0db1dd8b603a7e8d91e85e480d", label="ec_count")


tests = [test_mifaser_single_end_fastq]
