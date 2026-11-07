"""Tests for flyte_bio.modules.bedops.

The case and expected md5 are upstream's bedops/gtf2bed test (modules commit
9a1cc9f, BEDOPS 2.4.42).
"""

from flyte_bio.modules.bedops import bedops_gtf2bed
from tests.framework import assert_md5, env, fixture


@env.task
async def test_bedops_gtf2bed_homo_sapiens() -> None:
    bed = await bedops_gtf2bed(await fixture("genomics/homo_sapiens/genome/genome.gtf"))
    assert bed.name == "genome.bed", f"output is {bed.name}"
    await assert_md5(bed, "783257789b0c59a1a5391f7f2b9ac08a")


tests = [
    test_bedops_gtf2bed_homo_sapiens,
]
