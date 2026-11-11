"""Tests for flyte_bio.modules.wittyer.

The cases are upstream's wittyer tests (modules commit 9a1cc9f, witty.er
0.5.2): the decompressed simulated SV VCFs (query ``simulated_sv.vcf``, truth
``simulated_sv2.vcf``) compared plainly, within a BED and with a config. Each
checks the statistics JSON's md5, the index's name and nft-vcf's
``variantsMD5`` of the annotated VCF (the md5 of its record lines joined
without newlines). Upstream decompresses the fixtures with htslib's bgzip;
these tests use gunzip, which gives the same files.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.wittyer import wittyer
from tests.framework import assert_md5, env, fixture

SIMULATED_SV = "genomics/homo_sapiens/illumina/vcf/chr21/simulated_sv"
VARIANTS_MD5 = "de5a3d33e6ea3906502561b5f4d968f7"


async def wittyer_case(bed: File | None, config: File | None, report_md5: str) -> None:
    query, truth = await asyncio.gather(fixture(SIMULATED_SV + ".vcf.gz"), fixture(SIMULATED_SV + "2.vcf.gz"))
    query, truth = await asyncio.gather(gunzip(query), gunzip(truth))
    result = await wittyer(query, truth, prefix="test", bed=bed, config=config)
    await assert_md5(result.report, report_md5, label="report")
    assert result.bench_vcf_tbi.name == "test.vcf.gz.tbi", result.bench_vcf_tbi.name
    async with result.bench_vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    records = [x for x in lines if not x.startswith("#")]
    md5 = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
    assert md5 == VARIANTS_MD5, f"variantsMD5 {md5}"


@env.task
async def test_wittyer_vcf() -> None:
    await wittyer_case(None, None, "d1e6eccb85b0f4b5bf82336bbf1f8ce3")


@env.task
async def test_wittyer_vcf_bed() -> None:
    bed = await fixture("genomics/homo_sapiens/genome/chr21/sequence/multi_intervals.bed")
    await wittyer_case(bed, None, "8dfd8086fcb0ccc677fb31ea1f3dfb67")


@env.task
async def test_wittyer_vcf_config() -> None:
    config = await fixture("generic/config/wittyer_config.json")
    await wittyer_case(None, config, "c6361bea5cacdd463fe1a82966aca0bb")


tests = [test_wittyer_vcf, test_wittyer_vcf_bed, test_wittyer_vcf_config]
