"""Tests for flyte_bio.modules.pydamage.

The case and expected md5 are upstream's pydamage/analyze test (modules
commit 9a1cc9f, pydamage 1.0).
"""

from flyte_bio.modules.pydamage import pydamage_analyze
from tests.framework import assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"


@env.task
async def test_pydamage_analyze() -> None:
    # upstream case: "homo_sapiens"
    csv = await pydamage_analyze(await fixture(BAM), await fixture(BAM + ".bai"), prefix="test")
    await assert_md5(csv, "c0093b9583e99d3bac25cb35d8fe827d", label="test_pydamage_results.csv")


tests = [test_pydamage_analyze]
