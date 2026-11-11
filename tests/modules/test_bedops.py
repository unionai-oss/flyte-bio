"""Tests for flyte_bio.modules.bedops.

The case and expected md5 are upstream's bedops/gtf2bed test (modules commit
9a1cc9f, BEDOPS 2.4.42), and its bedops/convert2bed test at the same commit.
"""

from flyte_bio.modules.bedops import bedops_convert2bed, bedops_gtf2bed
from tests.framework import assert_md5, env, fixture


@env.task
async def test_bedops_gtf2bed_homo_sapiens() -> None:
    bed = await bedops_gtf2bed(await fixture("genomics/homo_sapiens/genome/genome.gtf"))
    assert bed.name == "genome.bed", f"output is {bed.name}"
    await assert_md5(bed, "783257789b0c59a1a5391f7f2b9ac08a")


@env.task
async def test_bedops_convert2bed_gtf() -> None:
    bed = await bedops_convert2bed(await fixture("genomics/sarscov2/genome/genome.gtf"), prefix="test")
    assert bed.name == "test.bed", f"output is {bed.name}"
    await assert_md5(bed, "adfdd36e848a62f4b0ea8a694abe9659")


tests = [
    test_bedops_gtf2bed_homo_sapiens,
    test_bedops_convert2bed_gtf,
]
