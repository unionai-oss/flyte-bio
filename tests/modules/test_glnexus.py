"""Tests for flyte_bio.modules.glnexus.

The cases and expected md5s are upstream's glnexus tests (modules commit
9a1cc9f, GLnexus 1.4.1). As upstream's module, the custom config is staged
but not passed, so all three cases write the same BCF.
"""

import asyncio

from flyte.io import File

from flyte_bio.modules.glnexus import glnexus
from tests.framework import assert_md5, env, fixture

GVCF = "genomics/homo_sapiens/illumina/gvcf/"
BCF_MD5 = "62b2cea9c1b92ac63645cb031eea46fc"


async def glnexus_gvcfs() -> list[File]:
    return list(await asyncio.gather(fixture(GVCF + "test.genome.vcf.gz"), fixture(GVCF + "test2.genome.vcf.gz")))


@env.task
async def test_glnexus_gvcfs() -> None:
    # upstream case: glnexus "vcfs, []"
    bcf = await glnexus(await glnexus_gvcfs(), prefix="test")
    assert bcf.path.endswith("/test.bcf"), bcf.path
    await assert_md5(bcf, BCF_MD5, label="test.bcf")


@env.task
async def test_glnexus_custom_config() -> None:
    # upstream case: glnexus "vcfs, [], custom_config"
    gvcfs, config = await asyncio.gather(glnexus_gvcfs(), fixture("delete_me/glnexus/gatk.yml"))
    bcf = await glnexus(gvcfs, prefix="test", custom_config=config)
    await assert_md5(bcf, BCF_MD5, label="test.bcf")


@env.task
async def test_glnexus_bed() -> None:
    # upstream case: glnexus "vcfs, bed"
    gvcfs, bed = await asyncio.gather(glnexus_gvcfs(), fixture("genomics/homo_sapiens/genome/genome.bed"))
    bcf = await glnexus(gvcfs, prefix="test", bed=bed)
    await assert_md5(bcf, BCF_MD5, label="test.bcf")


tests = [test_glnexus_gvcfs, test_glnexus_custom_config, test_glnexus_bed]
