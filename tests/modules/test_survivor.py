"""Tests for flyte_bio.modules.survivor.

The cases are upstream's survivor/filter, survivor/merge and survivor/stats
tests (modules commit 9a1cc9f, SURVIVOR 1.0.7): the filtered VCF's md5 (with
and without a BED of regions to ignore), nft-vcf's summary of the merged VCF
(contigs, samples, records) and the stats table's md5.
"""

import asyncio

from flyte.io import File

from flyte_bio.modules.survivor import survivor_filter, survivor_merge, survivor_stats
from tests.framework import assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
GVCF = HS + "illumina/gvcf/"


async def survivor_filter_case(bed: File | None) -> File:
    vcf = await fixture(GVCF + "test.genome.vcf")
    return await survivor_filter(vcf, "test", 51, 10001, 0.01, 10, bed=bed)


@env.task
async def test_survivor_filter_vcf() -> None:
    await assert_md5(await survivor_filter_case(None), "ec23f83ba5bb098d45027e63fab93d05")


@env.task
async def test_survivor_filter_vcf_bed() -> None:
    bed = await fixture(HS + "genome/genome.bed")
    await assert_md5(await survivor_filter_case(bed), "ec23f83ba5bb098d45027e63fab93d05")


@env.task
async def test_survivor_merge() -> None:
    vcfs = list(await asyncio.gather(fixture(GVCF + "test2.genome.vcf"), fixture(GVCF + "test.genome.vcf")))
    merged = await survivor_merge(vcfs, "test", 0.2, 1, 0, 0, 0, 0)
    assert merged.name == "test.vcf", merged.name
    async with merged.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    samples = len(next(x for x in lines if x.startswith("#CHROM")).split("\t")) - 9
    records = [x.split("\t") for x in lines if not x.startswith("#")]
    contigs = sorted({r[0] for r in records})
    summary = (contigs, samples, len(records))
    assert summary == (["chr22"], 2, 13), f"merged VCF (contigs, samples, records) {summary}"


@env.task
async def test_survivor_stats_vcfgz() -> None:
    vcf = await fixture(HS + "illumina/gatk/haplotypecaller_calls/test2_haplotc.vcf.gz")
    stats = await survivor_stats(vcf, "test", -1, -1, -1)
    await assert_md5(stats, "d51fb5cc28255f37096cd0f3a6de5f8d")


tests = [test_survivor_filter_vcf, test_survivor_filter_vcf_bed, test_survivor_merge, test_survivor_stats_vcfgz]
