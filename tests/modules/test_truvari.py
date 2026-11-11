"""Tests for flyte_bio.modules.truvari.

The cases are upstream's truvari/bench tests (modules commit 9a1cc9f,
Truvari 5.4.0), with and without a BED of regions: nft-vcf's summary of each
output VCF (contigs, samples, records), the indexes' and log's names, and the
summary JSON's md5.
"""

import asyncio
import gzip

from flyte.io import File

from flyte_bio.modules.truvari import TruvariBenchResult, truvari_bench
from tests.framework import assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
SIMULATED_SV = HS + "illumina/vcf/chr21/simulated_sv"
CHR21 = HS + "genome/chr21/sequence/genome.fasta"


async def truvari_vcf_summary(vcf: File) -> tuple[list[str], int, int]:
    """A gzipped VCF's contigs (in order of first record), sample count and record count."""
    async with vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    samples = len(next(x for x in lines if x.startswith("#CHROM")).split("\t")) - 9
    records = [x.split("\t", 1)[0] for x in lines if not x.startswith("#")]
    return list(dict.fromkeys(records)), samples, len(records)


async def truvari_case(bed: File | None) -> TruvariBenchResult:
    vcf, tbi, truth, truth_tbi, fasta, fai = await asyncio.gather(
        fixture(SIMULATED_SV + ".vcf.gz"),
        fixture(SIMULATED_SV + ".vcf.gz.tbi"),
        fixture(SIMULATED_SV + "2.vcf.gz"),
        fixture(SIMULATED_SV + "2.vcf.gz.tbi"),
        fixture(CHR21),
        fixture(CHR21 + ".fai"),
    )
    return await truvari_bench(vcf, tbi, truth, truth_tbi, fasta, fai, prefix="test", bed=bed)


async def check(result: TruvariBenchResult, records: dict[str, int], summary_md5: str) -> None:
    for name, count in records.items():
        key = name.replace("-", "_")
        vcf, tbi = getattr(result, f"{key}_vcf"), getattr(result, f"{key}_tbi")
        assert tbi.name == f"test.{name}.vcf.gz.tbi", tbi.name
        expected = (["chr21"] if count else [], 1, count)
        actual = await truvari_vcf_summary(vcf)
        assert actual == expected, f"{name}: (contigs, samples, records) {actual}, expected {expected}"
    assert result.log.name == "test.log.txt", result.log.name
    await assert_md5(result.summary, summary_md5, label="summary")


@env.task
async def test_truvari_bench_without_bed() -> None:
    result = await truvari_case(None)
    await check(result, {"fn": 136, "fp": 136, "tp-base": 0, "tp-comp": 0}, "aff467e3dbccf2c00a7aa81c2a36f62a")


@env.task
async def test_truvari_bench_with_bed() -> None:
    result = await truvari_case(await fixture(HS + "genome/genome.bed"))
    await check(result, {"fn": 0, "fp": 0, "tp-base": 0, "tp-comp": 0}, "397452f6fd093809aee763a3d2310552")


tests = [test_truvari_bench_without_bed, test_truvari_bench_with_bed]
