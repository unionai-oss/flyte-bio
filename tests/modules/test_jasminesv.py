"""Tests for flyte_bio.modules.jasminesv.

The cases are upstream's jasminesv tests (modules commit 9a1cc9f, Jasmine
1.1.5), with their params' ``args`` and ``args2`` (Iris options). As
upstream, each checks the merged VCF's summary (chromosomes, samples,
variant count, phasing); Jasmine's header records the merge command.
"""

import asyncio
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.jasminesv import jasminesv
from tests.framework import env, fixture
from tests.modules.test_cutesv import VcfSummary, vcf_summary

HUMAN = "genomics/homo_sapiens/"
GVCF = HUMAN + "illumina/gvcf/"
BAM = HUMAN + "illumina/bam/"
IRIS_ARGS = {"args": "--run_iris", "args2": "threads=1,--pacbio,min_ins_length=30"}


async def jasminesv_gvcfs() -> list[File]:
    return list(await asyncio.gather(fixture(GVCF + "test.genome.vcf"), fixture(GVCF + "test2.genome.vcf")))


async def jasminesv_bams() -> list[File]:
    return list(
        await asyncio.gather(fixture(BAM + "test.paired_end.sorted.bam"), fixture(BAM + "test2.paired_end.sorted.bam"))
    )


async def jasminesv_reference() -> tuple[File, File]:
    fasta, fai = await asyncio.gather(
        fixture(HUMAN + "genome/genome.fasta"), fixture(HUMAN + "genome/genome.fasta.fai")
    )
    return fasta, fai


async def check(vcf: File, expected: VcfSummary) -> None:
    assert vcf.path.endswith("/test.vcf.gz"), vcf.path
    summary = await vcf_summary(vcf)
    assert summary == expected, summary


@env.task
async def test_jasminesv_vcf() -> None:
    # upstream case: jasminesv "homo_sapiens - [vcf, vcf2], [], [], [], [], [], [] - vcf"
    vcf = await jasminesv(await jasminesv_gvcfs(), prefix="test")
    await check(vcf, VcfSummary(["chr22"], 1, 207, False))


@env.task
async def test_jasminesv_vcf_gz() -> None:
    # upstream case: jasminesv "homo_sapiens - [vcf_gz, vcf_gz2], [], [], [], [], [], [] - vcf.gz"
    calls = HUMAN + "illumina/gatk/haplotypecaller_calls/"
    vcfs = list(
        await asyncio.gather(
            fixture(calls + "test2_haplotc.ann.vcf.gz"), fixture(calls + "test_haplotcaller.cnn.vcf.gz")
        )
    )
    vcf = await jasminesv(vcfs, prefix="test")
    await check(vcf, VcfSummary(["chr21", "chr22"], 1, 738, False))


@env.task
async def test_jasminesv_iris() -> None:
    # upstream case: jasminesv "homo_sapiens - [vcf, vcf2], [bam, bam2], [], [], fasta, fai, [] - iris"
    vcfs, bams, (fasta, fai) = await asyncio.gather(jasminesv_gvcfs(), jasminesv_bams(), jasminesv_reference())
    vcf = await jasminesv(vcfs, prefix="test", bams=bams, fasta=fasta, fai=fai, **IRIS_ARGS)
    await check(vcf, VcfSummary(["chr22"], 1, 205, False))


@env.task
async def test_jasminesv_iris_indexed() -> None:
    # upstream case: jasminesv "homo_sapiens - [vcf, vcf2], [bam, bam2], [bai1, bai2], [], fasta, fai, [] - iris"
    vcfs, bams, (fasta, fai), bais = await asyncio.gather(
        jasminesv_gvcfs(),
        jasminesv_bams(),
        jasminesv_reference(),
        asyncio.gather(
            fixture(BAM + "test.paired_end.sorted.bam.bai"), fixture(BAM + "test2.paired_end.sorted.bam.bai")
        ),
    )
    vcf = await jasminesv(vcfs, prefix="test", bams=bams, bais=list(bais), fasta=fasta, fai=fai, **IRIS_ARGS)
    await check(vcf, VcfSummary(["chr22"], 1, 205, False))


@env.task
async def test_jasminesv_normalize() -> None:
    # upstream case: jasminesv "homo_sapiens - [vcf, vcf2], [bam, bam2], [], fasta, fai, txt - normalize"
    path = Path(tempfile.mkdtemp()) / "chr_norm.txt"
    path.write_text("chr21 21\nchr22 22\n")
    vcfs, bams, (fasta, fai), chr_norm = await asyncio.gather(
        jasminesv_gvcfs(), jasminesv_bams(), jasminesv_reference(), File.from_local(str(path))
    )
    vcf = await jasminesv(
        vcfs, prefix="test", bams=bams, fasta=fasta, fai=fai, chr_norm=chr_norm, args="--normalize_chrs"
    )
    await check(vcf, VcfSummary(["22"], 1, 207, False))


@env.task
async def test_jasminesv_empty() -> None:
    # upstream case: jasminesv "... - empty vcf - not stub, but linting complains otherwise"
    vcfs = list(await asyncio.gather(fixture(GVCF + "test.genome.vcf"), fixture(HUMAN + "illumina/vcf/empty.vcf.gz")))
    vcf = await jasminesv(vcfs, prefix="test", args="min_support=2")
    await check(vcf, VcfSummary([], 1, 0, True))


tests = [
    test_jasminesv_vcf,
    test_jasminesv_vcf_gz,
    test_jasminesv_iris,
    test_jasminesv_iris_indexed,
    test_jasminesv_normalize,
    test_jasminesv_empty,
]
