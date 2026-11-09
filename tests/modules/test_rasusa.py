"""Tests for flyte_bio.modules.rasusa.

The case and expected md5s are upstream's rasusa test (modules commit
9a1cc9f, Rasusa 0.3.0): the paired SARS-CoV-2 reads at 100x of a 1 Mb
genome. Gzipped outputs are compared decompressed, as nf-test hashes them.
"""

from flyte_bio.modules.rasusa import rasusa
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_rasusa_paired_end() -> None:
    # upstream case: "Should run without failures"
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    r1, r2 = await rasusa(reads, "1000000b", 100, prefix="testfile")
    await assert_gunzipped_md5(r1, "4161df271f9bfcd25d5845a1e220dbec", label="testfile_1.fastq.gz")
    await assert_gunzipped_md5(r2, "2ebae722295ea66d84075a3b042e2b42", label="testfile_2.fastq.gz")


tests = [test_rasusa_paired_end]
