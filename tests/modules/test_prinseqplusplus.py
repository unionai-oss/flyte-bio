"""Tests for flyte_bio.modules.prinseqplusplus.

The cases are upstream's prinseqplusplus tests (modules commit 9a1cc9f,
PRINSEQ++ 1.2.3), with the test directory's ``-lc_entropy=0.8``. Upstream
marks every output unstable and snapshots their names only, as these tests
check (plus that the passing reads aren't empty).
"""

from flyte.io import File

from flyte_bio.modules.prinseqplusplus import prinseqplusplus
from tests.framework import assert_nonempty, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
ARGS = "-lc_entropy=0.8"


def names(files: list[File]) -> list[str]:
    return [f.name for f in files]


@env.task
async def test_prinseqplusplus_single_end() -> None:
    # upstream case: prinseqplusplus "test-prinseqplusplus-single-end"
    result = await prinseqplusplus(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", args=ARGS)
    assert names(result.good_reads) == ["test_good_out.fastq.gz"], names(result.good_reads)
    assert names(result.bad_reads) == ["test_bad_out.fastq.gz"], names(result.bad_reads)
    assert result.single_reads == [], names(result.single_reads)
    assert result.log.name == "test.log", result.log.name
    await assert_nonempty(result.good_reads[0], label="test_good_out.fastq.gz")


@env.task
async def test_prinseqplusplus_paired_end() -> None:
    # upstream case: prinseqplusplus "test-prinseqplusplus-paired-end"
    result = await prinseqplusplus(
        await fixture(FASTQ + "test_1.fastq.gz"),
        reads_2=await fixture(FASTQ + "test_2.fastq.gz"),
        prefix="test",
        args=ARGS,
    )
    assert names(result.good_reads) == ["test_good_out_R1.fastq.gz", "test_good_out_R2.fastq.gz"]
    assert names(result.single_reads) == ["test_single_out_R1.fastq.gz", "test_single_out_R2.fastq.gz"]
    assert names(result.bad_reads) == ["test_bad_out_R1.fastq.gz", "test_bad_out_R2.fastq.gz"]
    assert result.log.name == "test.log", result.log.name
    for read in result.good_reads:
        await assert_nonempty(read, label=read.name)


tests = [test_prinseqplusplus_single_end, test_prinseqplusplus_paired_end]
