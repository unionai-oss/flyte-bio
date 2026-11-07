"""Tests for flyte_bio.modules.clair3.

The cases and expected values are upstream's clair3 tests (modules commit
9a1cc9f, Clair3 2.0.0), all with ``--include_all_ctgs`` and the HiFi
platform: the SARS-CoV-2 test BAM with a user model (upstream's
``clair3_v2.hifi.tar.gz``, extracted, as its test's untar step), with the
image's packaged ``hifi_revio`` model, with both (which must fail), and with a
user model and ``--gvcf``. Upstream snapshots each VCF's nft-vcf summary
(contigs, sample count, record count, phasing) and the indexes' names; these
tests check the same.
"""

import asyncio
import gzip

from flyte.io import Dir, File

from flyte_bio.modules.clair3 import Clair3Result, clair3
from tests.framework import env, fixture, fixture_dir

BAM = "genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"
GENOME = "genomics/sarscov2/genome/genome.fasta"
MODEL = "generic/models/clair3_v2.hifi.tar.gz"
ARGS = "--include_all_ctgs"


async def vcf_summary(vcf: File) -> dict:
    """What nft-vcf's ``summary`` reports: contigs, samples, records and whether any genotype is phased."""
    local = await vcf.download()
    chromosomes: list[str] = []
    samples = records = 0
    phased = False
    with gzip.open(local, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = len(fields) - 9
                continue
            records += 1
            if fields[0] not in chromosomes:
                chromosomes.append(fields[0])
            if len(fields) > 9 and "|" in fields[9].split(":")[0]:
                phased = True
    return {"chromosomes": chromosomes, "samples": samples, "records": records, "phased": phased}


async def check(result: Clair3Result, gvcf_records: int | None = None) -> None:
    assert result.vcf.name == "testmerge_output.vcf.gz", result.vcf.name
    assert result.tbi.name == "testmerge_output.vcf.gz.tbi", result.tbi.name
    expected = {"chromosomes": ["MT192765.1"], "samples": 1, "records": 9, "phased": False}
    summary = await vcf_summary(result.vcf)
    assert summary == expected, f"vcf summary {summary} != {expected}"
    assert result.phased_vcf is None and result.phased_tbi is None, result
    if gvcf_records is None:
        assert result.gvcf is None and result.gvcf_tbi is None, result
    else:
        assert result.gvcf is not None and result.gvcf.name == "testmerge_output.gvcf.gz", result.gvcf
        assert result.gvcf_tbi is not None and result.gvcf_tbi.name == "testmerge_output.gvcf.gz.tbi", result
        summary = await vcf_summary(result.gvcf)
        expected = {**expected, "records": gvcf_records}
        assert summary == expected, f"gvcf summary {summary} != {expected}"


async def inputs() -> tuple[File, File, File, File]:
    return await asyncio.gather(fixture(BAM), fixture(BAM + ".bai"), fixture(GENOME), fixture(GENOME + ".fai"))


async def user_model() -> Dir:
    return await fixture_dir(MODEL, "hifi_revio")


@env.task
async def test_clair3_user_model() -> None:
    # upstream case: clair3 "sarscov2 - bam - user_model"
    (bam, bai, fasta, fai), model = await asyncio.gather(inputs(), user_model())
    result = await clair3(bam, bai, fasta, fai, "hifi", user_model=model, prefix="test", args=ARGS)
    await check(result)


@env.task
async def test_clair3_packaged_model() -> None:
    # upstream case: clair3 "sarscov2 - bam - packaged_model"
    bam, bai, fasta, fai = await inputs()
    result = await clair3(bam, bai, fasta, fai, "hifi", packaged_model="hifi_revio", prefix="test", args=ARGS)
    await check(result)


@env.task
async def test_clair3_both_models() -> None:
    # upstream case: clair3 "sarscov2 - bam - both" must fail: give one model.
    (bam, bai, fasta, fai), model = await asyncio.gather(inputs(), user_model())
    try:
        await clair3(bam, bai, fasta, fai, "hifi", packaged_model="hifi_revio", user_model=model, prefix="test")
    except ValueError as e:
        assert "specify one of them" in str(e), e
    else:
        raise AssertionError("clair3 accepted two models")


@env.task
async def test_clair3_user_model_gvcf() -> None:
    # upstream case: clair3 "sarscov2 - bam - user_model - gvcf"
    (bam, bai, fasta, fai), model = await asyncio.gather(inputs(), user_model())
    result = await clair3(bam, bai, fasta, fai, "hifi", user_model=model, prefix="test", args="--gvcf " + ARGS)
    await check(result, gvcf_records=456)


tests = [test_clair3_user_model, test_clair3_packaged_model, test_clair3_both_models, test_clair3_user_model_gvcf]
