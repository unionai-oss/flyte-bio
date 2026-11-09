"""Tests for flyte_bio.modules.vt.

The cases are upstream's vt/decompose and vt/normalize tests (modules commit
9a1cc9f, vt 0.57721): decompose's decompressed output md5, and normalize's
calls by nft-vcf's ``variantsMD5`` (or, within intervals, its summary: no
records, one sample) with the reference's ``.fai`` by md5.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.vt import VtNormalizeResult, vt_decompose, vt_normalize
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

SC2 = "genomics/sarscov2/"
HS = "genomics/homo_sapiens/"
CALLS = HS + "illumina/gatk/haplotypecaller_calls/test_haplotcaller.cnn.vcf.gz"


@env.task
async def test_vt_decompose() -> None:
    # upstream cases: vt/decompose "test-vt-decompose" and "test-vt-decompose-intervals"
    vcf_gz, vcf, bed = await asyncio.gather(
        fixture(SC2 + "illumina/vcf/test.vcf.gz"),
        fixture(SC2 + "illumina/vcf/test.vcf"),
        fixture(SC2 + "genome/bed/test.bed"),
    )
    plain, within = await asyncio.gather(
        vt_decompose(vcf_gz, prefix="test_decompose"), vt_decompose(vcf, prefix="test_decompose", intervals=bed)
    )
    for out, label in ((plain, "decompose"), (within, "decompose - intervals")):
        assert out.name == "test_decompose.vcf.gz", f"{label}: {out.name}"
        await assert_gunzipped_md5(out, "892b6b5c55c69afd3a96cd9bfd8ab743", label=label)


async def vcf_lines(vcf: File) -> tuple[list[str], list[str]]:
    """A gzipped VCF's (header lines, record lines)."""
    async with vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    return [x for x in lines if x.startswith("#")], [x for x in lines if not x.startswith("#")]


async def check_normalized(result: VtNormalizeResult, label: str) -> list[str]:
    assert result.vcf.name == "test_normalize.vcf.gz", f"{label}: {result.vcf.name}"
    assert result.fai is not None and result.fai.name == "genome.fasta.fai", f"{label}: {result.fai}"
    await assert_md5(result.fai, "3520cd30e1b100e55f578db9c855f685", label=f"{label} fai")
    header, records = await vcf_lines(result.vcf)
    columns = header[-1].split("\t")
    assert len(columns) == 10, f"{label}: {len(columns) - 9} samples, expected 1"
    return records


@env.task
async def test_vt_normalize() -> None:
    # upstream cases: vt/normalize "test-vt-normalize", "test-vt-normalize-intervals" and "test-vt-normalize-no-fai"
    vcf, tbi, fasta, fai, bed = await asyncio.gather(
        fixture(CALLS),
        fixture(CALLS + ".tbi"),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/genome.fasta.fai"),
        fixture(HS + "genome/genome.bed"),
    )
    common = {"vcf": vcf, "tbi": tbi, "fasta": fasta, "prefix": "test_normalize"}
    with_fai, within, without_fai = await asyncio.gather(
        vt_normalize(fai=fai, **common),
        vt_normalize(fai=fai, intervals=bed, **common),
        vt_normalize(**common),
    )
    for result, label in ((with_fai, "normalize"), (without_fai, "normalize - no fai")):
        records = await check_normalized(result, label)
        md5 = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
        assert md5 == "cdd99bce7ca89d001a65d448be3b4bae", f"{label}: variantsMD5 {md5}"
    records = await check_normalized(within, "normalize - intervals")
    assert records == [], f"normalize - intervals: {len(records)} records, expected none"


tests = [test_vt_decompose, test_vt_normalize]
