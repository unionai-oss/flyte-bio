"""Tests for flyte_bio.modules.subread.

Expected md5s are the upstream subread/featurecounts snapshot values. The
counts header records the command line, including ``-T <threads>``; upstream's
snapshots were made with 2 threads (verified: ``-T 2`` reproduces them
exactly, other counts don't), so these pass ``threads=2``.
"""

from flyte_bio.modules.subread import featurecounts
from tests.framework import assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/test.single_end.bam"
GTF = "genomics/sarscov2/genome/genome.gtf"


async def check(strandedness: str, counts_md5: str, summary_md5: str) -> None:
    bam, gtf = await fixture(BAM), await fixture(GTF)
    result = await featurecounts(
        bam, gtf, prefix="test", strandedness=strandedness, single_end=True, threads=2, args="-t CDS"
    )
    await assert_md5(result.counts, counts_md5, label=f"{strandedness} counts")
    await assert_md5(result.summary, summary_md5, label=f"{strandedness} summary")


@env.task
async def test_featurecounts_forward() -> None:
    # upstream case: subread/featurecounts "sarscov2 [bam] - forward"
    await check("forward", "21ff44bfaa4a3d8e8b7e749078f7a201", "8f602ff9a8ef467af43294e80b367cdf")


@env.task
async def test_featurecounts_reverse() -> None:
    # upstream case: subread/featurecounts "sarscov2 [bam] - reverse"
    await check("reverse", "7df6092fdc65ce40b71c64e0c97f95c6", "7cfa30ad678b9bc1bc63afbb0281547b")


@env.task
async def test_featurecounts_unstranded() -> None:
    # upstream case: subread/featurecounts "sarscov2 [bam] - unstranded"
    await check("unstranded", "9474f78d2d1d43613cbc16c10ba15047", "23164b79f9f23f11c82820db61a35560")


tests = [test_featurecounts_forward, test_featurecounts_reverse, test_featurecounts_unstranded]
