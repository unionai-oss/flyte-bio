"""Tests for flyte_bio.modules.seacr.

The cases and expected md5s are upstream's seacr/callpeak tests (modules
commit 9a1cc9f, SEACR 1.3), run with upstream's test ``args``,
``norm stringent``.
"""

from flyte_bio.modules.seacr import seacr_callpeak
from tests.framework import assert_md5, env, fixture

BEDGRAPH = "genomics/homo_sapiens/illumina/bedgraph/"
ARGS = "norm stringent"


@env.task
async def test_seacr_callpeak_control() -> None:
    target = await fixture(BEDGRAPH + "cutandtag_h3k27me3_test_1.bedGraph")
    control = await fixture(BEDGRAPH + "cutandtag_igg_test_1.bedGraph")
    bed = await seacr_callpeak(target, control, threshold=0.05, prefix="test_1", args=ARGS)
    assert bed.name == "test_1.stringent.bed", f"seacr: wrote {bed.name}"
    await assert_md5(bed, "a3cb0c7c4ffa895788da3f0d6371b7df", label=bed.name)


@env.task
async def test_seacr_callpeak_threshold() -> None:
    target = await fixture(BEDGRAPH + "cutandtag_h3k27me3_test_1.bedGraph")
    bed = await seacr_callpeak(target, None, threshold=0.05, prefix="test_1", args=ARGS)
    assert bed.name == "test_1.stringent.bed", f"seacr: wrote {bed.name}"
    await assert_md5(bed, "1d23015c7087f7b48cc3139d53fd3463", label=bed.name)


tests = [test_seacr_callpeak_control, test_seacr_callpeak_threshold]
