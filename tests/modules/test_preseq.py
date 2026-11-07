"""Tests for flyte_bio.modules.preseq.

The cases and expected md5s are upstream's preseq/lcextrap tests (modules
commit 9a1cc9f, preseq 3.2.0). preseq's bootstrap uses a fixed default seed,
so the extrapolation is reproducible. As upstream, the log is only checked by
name (preseq writes nothing to stderr without ``-v``, so it may be empty).

The c_curve case is upstream's preseq/ccurve tests (same commit). Upstream's
"single-end" and "paired-end" cases are identical (both run the paired-end
BAM with ``single_end: false``) and share a snapshot, so one test covers both.
"""

from flyte_bio.modules.preseq import preseq_ccurve, preseq_lcextrap
from tests.framework import assert_md5, env, fixture

MAPPED_READS = "delete_me/preseq/SRR1003759_5M_subset.mr"


@env.task
async def test_preseq_lcextrap_single_end() -> None:
    result = await preseq_lcextrap(await fixture(MAPPED_READS), prefix="test", single_end=True)
    await assert_md5(result.lc_extrap, "1fa5cdd601079329618f61660bee00de", label="test.lc_extrap.txt")
    assert result.log.name == "test.command.log", f"preseq log: {result.log.name}"


@env.task
async def test_preseq_lcextrap_paired_end() -> None:
    result = await preseq_lcextrap(await fixture(MAPPED_READS), prefix="test", single_end=False)
    await assert_md5(result.lc_extrap, "10e5ea860e87fb6f5dc10f4f20c62040", label="test.lc_extrap.txt")
    assert result.log.name == "test.command.log", f"preseq log: {result.log.name}"


@env.task
async def test_preseq_ccurve_paired_end() -> None:
    # upstream cases: preseq/ccurve "test-preseq-ccurve-single-end", "test-preseq-ccurve-paired-end"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.bam")
    result = await preseq_ccurve(bam, prefix="test", single_end=False)
    await assert_md5(result.c_curve, "cf4743abdd355595d6ec1fb3f38e66e5", label="test.c_curve.txt")
    assert result.log.name == "test.command.log", f"preseq log: {result.log.name}"


tests = [test_preseq_lcextrap_single_end, test_preseq_lcextrap_paired_end, test_preseq_ccurve_paired_end]
