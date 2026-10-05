"""Tests for flyte_bio.modules.bbmap.

Mirrors the upstream bbmap/bbsplit "split with prebuilt index" case: build
an index, then split reads with it. Expected md5s are the upstream snapshot
values (nf-test hashes the decompressed reads). The log check proves the
staged index was reused rather than rebuilt — the timestamp/path fix-up
upstream applies to a staged index is what makes that work.
"""

from flyte_bio.modules.bbmap import bbsplit, bbsplit_index
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture


@env.task
async def test_bbsplit_prebuilt_index() -> None:
    # upstream case: bbmap/bbsplit "sarscov2_se_fastq_fasta_chr22_fasta - split with prebuilt index"
    primary = await fixture("genomics/sarscov2/genome/genome.fasta")
    human = await fixture("genomics/homo_sapiens/genome/chr22/sequence/chr22_23800000-23980000.fa")
    index = await bbsplit_index(primary, {"human": human})

    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    result = await bbsplit(reads, None, index, prefix="test")
    await assert_gunzipped_md5(result.reads_1, "4161df271f9bfcd25d5845a1e220dbec", label="primary reads")
    await assert_md5(result.stats, "2cbf69b72e5f4f8508306b54e8fe2861", label="stats")
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    assert "If you wish to regenerate the index" in log, "bbsplit rebuilt the index instead of reusing it"


tests = [test_bbsplit_prebuilt_index]
