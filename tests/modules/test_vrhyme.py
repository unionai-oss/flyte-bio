"""Tests for flyte_bio.modules.vrhyme.

The case is upstream's vrhyme/vrhyme test (modules commit 9a1cc9f, vRhyme
1.1.0): the Bacteroides fragilis contigs binned with its paired-end reads.
Upstream snapshots the bins' names and the membership table's md5, and
the summary's name only, as this test checks. The module's 4 CPUs match
upstream's test CPUs.
"""

from flyte_bio.modules.vrhyme import vrhyme
from tests.framework import assert_md5, env, fixture

BACTEROIDES = "genomics/prokaryotes/bacteroides_fragilis/illumina/"
BIN_NAMES = sorted(f"vRhyme_bin_{n}.{ext}" for n in range(1, 15) for ext in ("faa", "fasta", "ffn"))


@env.task
async def test_vrhyme_contigs() -> None:
    # upstream case: vrhyme/vrhyme "test-vrhyme-vrhyme"
    result = await vrhyme(
        await fixture(BACTEROIDES + "fastq/test1_1.fastq.gz"),
        await fixture(BACTEROIDES + "fasta/test1.contigs.fa.gz"),
        reads_2=await fixture(BACTEROIDES + "fastq/test1_2.fastq.gz"),
        prefix="test",
    )
    names = sorted(f.name for f in await result.bins.list_files())
    assert names == BIN_NAMES, f"bins: {names}"
    assert result.membership.name == "vRhyme_best_bins.19.membership.tsv", result.membership.name
    await assert_md5(result.membership, "a24fe8417dd21c2db3da84cf51cff750", label=result.membership.name)
    assert result.summary.name == "vRhyme_best_bins.19.summary.tsv", result.summary.name


tests = [test_vrhyme_contigs]
