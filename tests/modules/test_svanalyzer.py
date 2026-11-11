"""Tests for flyte_bio.modules.svanalyzer.

The cases are upstream's svanalyzer/svbenchmark tests (modules commit
9a1cc9f, SVanalyzer 0.36), with and without a BED of regions: the md5s of the
false-negative and false-positive VCFs (decompressed, as nf-test hashes them),
the distances and the report. (Upstream's snapshot names are swapped: the case
with a BED is stored under ``fns``…, the one without under ``bed_fns``….)
"""

import asyncio

from flyte.io import File

from flyte_bio.modules.svanalyzer import SvbenchmarkResult, svanalyzer_svbenchmark
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
CALLS = HS + "illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz"
TRUTH = HS + "illumina/gatk/haplotypecaller_calls/test2_haplotc.vcf.gz"
CHR21 = HS + "genome/chr21/sequence/"


async def svanalyzer_case(bed: File | None) -> SvbenchmarkResult:
    vcf, tbi, truth, truth_tbi, fasta, fai = await asyncio.gather(
        fixture(CALLS),
        fixture(CALLS + ".tbi"),
        fixture(TRUTH),
        fixture(TRUTH + ".tbi"),
        fixture(CHR21 + "genome.fasta"),
        fixture(CHR21 + "genome.fasta.fai"),
    )
    return await svanalyzer_svbenchmark(vcf, tbi, truth, truth_tbi, fasta, fai, prefix="test", bed=bed)


async def check(result: SvbenchmarkResult, fns: str, fps: str, distances: str, report: str) -> None:
    await assert_gunzipped_md5(result.fns, fns, label="false negatives")
    await assert_gunzipped_md5(result.fps, fps, label="false positives")
    await assert_md5(result.distances, distances, label="distances")
    await assert_md5(result.report, report, label="report")


@env.task
async def test_svanalyzer_svbenchmark_bed() -> None:
    result = await svanalyzer_case(await fixture(CHR21 + "multi_intervals.bed"))
    await check(
        result,
        fns="c04105833d5f9acb5ec50af0b89a6b0f",
        fps="146f05125330a9c0c4581a5a581744a7",
        distances="95f7b1679e805979f0c9a4a322cb77f4",
        report="aaa9225dd6e261951c020e23473375b9",
    )


@env.task
async def test_svanalyzer_svbenchmark_nobed() -> None:
    result = await svanalyzer_case(None)
    await check(
        result,
        fns="7a574ad3c614f0c38a6dccb5901ddc05",
        fps="146f05125330a9c0c4581a5a581744a7",
        distances="b6a908b1187f159429dc8c6beb6e6f57",
        report="d0c7db2466442cbdd4f87ad516e00d32",
    )


tests = [test_svanalyzer_svbenchmark_bed, test_svanalyzer_svbenchmark_nobed]
