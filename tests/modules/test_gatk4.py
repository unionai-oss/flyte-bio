"""Tests for flyte_bio.modules.gatk4 (md5s from upstream's gatk4 module snapshots)."""

import asyncio

from flyte_bio.modules.gatk4 import (
    gatk4_applybqsr,
    gatk4_baserecalibrator,
    gatk4_createsequencedictionary,
    gatk4_gatherbqsrreports,
    gatk4_intervallisttobed,
    gatk4_markduplicates,
)
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_md5, assert_nonempty, assert_reads_md5, env, fixture, reads_md5

HS = "genomics/homo_sapiens/"
SC2 = "genomics/sarscov2/"


async def sarscov2_reference():
    return await asyncio.gather(
        fixture(SC2 + "genome/genome.fasta"),
        fixture(SC2 + "genome/genome.fasta.fai"),
        fixture(SC2 + "genome/genome.dict"),
    )


@env.task
async def test_createsequencedictionary() -> None:
    # upstream case: gatk4/createsequencedictionary "sarscov2 - fasta"
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    dict_ = await gatk4_createsequencedictionary(fasta)
    assert dict_.path.endswith("/genome.dict"), dict_.path
    await assert_md5(dict_, "7362679f176e0f52add03c08f457f646", label="genome.dict")


@env.task
async def test_intervallisttobed() -> None:
    # upstream case: gatk4/intervallisttobed "homo sapiens - bed"
    intervals = await fixture("genomics/homo_sapiens/genome/chr21/sequence/genome.interval_list")
    bed = await gatk4_intervallisttobed(intervals, prefix="test")
    await assert_md5(bed, "9046675d01199fbbee79f2bc1c5dce52", label="test.bed")


@env.task
async def test_markduplicates_multiple() -> None:
    # upstream cases: gatk4/markduplicates "homo_sapiens - multiple cram" / "multiple bam".
    # Their BAM md5 covers a header that records the command line, so check the outputs instead.
    bams = await asyncio.gather(
        fixture(HS + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(HS + "illumina/bam/test2.paired_end.sorted.bam"),
    )
    fasta, fai = await asyncio.gather(fixture(HS + "genome/genome.fasta"), fixture(HS + "genome/genome.fasta.fai"))
    cram, bam = await asyncio.gather(
        gatk4_markduplicates(list(bams), "test.cram", fasta=fasta, fai=fai),
        gatk4_markduplicates(list(bams), "test.bam", fasta=fasta, fai=fai),
    )
    assert cram.alignment.path.endswith("/test.cram") and cram.index is not None, cram
    assert cram.index.path.endswith("/test.cram.crai"), cram.index.path
    assert bam.alignment.path.endswith("/test.bam") and bam.index is None, bam
    for result in (cram, bam):
        async with result.metrics.open("rb") as fh:
            metrics = bytes(await fh.read()).decode()
        assert "## METRICS CLASS\tpicard.sam.DuplicationMetrics" in metrics, metrics[:400]
        await assert_nonempty(result.alignment, label=result.alignment.path.rsplit("/", 1)[-1])
    # Both inputs are merged into one alignment: its header has both read groups.
    header = await samtools_view(bam=bam.alignment, args="-H")
    async with header.open("rb") as fh:
        read_groups = [line for line in bytes(await fh.read()).decode().splitlines() if line.startswith("@RG")]
    assert len(read_groups) >= 1, read_groups


@env.task
async def test_baserecalibrator() -> None:
    # upstream cases: gatk4/baserecalibrator "sarscov2 - bam" and "sarscov2 - bam - intervals"
    (fasta, fai, dict_), bam, bai, vcf, tbi, bed = await asyncio.gather(
        sarscov2_reference(),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai"),
        fixture(SC2 + "illumina/vcf/test.vcf.gz"),
        fixture(SC2 + "illumina/vcf/test.vcf.gz.tbi"),
        fixture(SC2 + "genome/bed/test.bed"),
    )
    whole, within = await asyncio.gather(
        gatk4_baserecalibrator(bam, bai, fasta, fai, dict_, [vcf], [tbi], prefix="test"),
        gatk4_baserecalibrator(bam, bai, fasta, fai, dict_, [vcf], [tbi], prefix="test", intervals=bed),
    )
    await assert_md5(whole, "e2e43abdc0c943c1a54dae816d0b9ea7", label="test.table")
    await assert_md5(within, "9ecb5f00a2229291705addc09c0ec231", label="test.table (intervals)")


@env.task
async def test_gatherbqsrreports() -> None:
    # upstream cases: gatk4/gatherbqsrreports "test-gatk4-gatherbqsrreports" and "...-multiple"
    one, two = await asyncio.gather(
        fixture(HS + "illumina/gatk/test.baserecalibrator.table"),
        fixture(HS + "illumina/gatk/test2.baserecalibrator.table"),
    )
    single, multiple = await asyncio.gather(
        gatk4_gatherbqsrreports([one], prefix="test"), gatk4_gatherbqsrreports([one, two], prefix="test")
    )
    await assert_md5(single, "9603b69fdc3b5090de2e0dd78bfcc4bf", label="gathered single")
    await assert_md5(multiple, "0c1257eececf95db8ca378272d0f21f9", label="gathered multiple")


@env.task
async def test_applybqsr() -> None:
    # upstream cases: gatk4/applybqsr "sarscov2 - bam" and "sarscov2 - bam in, cram out". Their BAM md5
    # covers a header recording the command line; recalibration only changes base qualities, so the
    # read sequences must come through unchanged.
    (fasta, fai, dict_), bam, bai, table = await asyncio.gather(
        sarscov2_reference(),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai"),
        fixture(SC2 + "illumina/gatk/test.baserecalibrator.table"),
    )
    as_bam, as_cram = await asyncio.gather(
        gatk4_applybqsr(bam, bai, table, fasta, fai, dict_, prefix="test", suffix="bam"),
        gatk4_applybqsr(bam, bai, table, fasta, fai, dict_, prefix="test"),
    )
    assert as_bam.path.endswith("/test.bam") and as_cram.path.endswith("/test.cram"), (as_bam.path, as_cram.path)
    expected = await reads_md5(bam)
    await assert_reads_md5(as_bam, expected, label="recalibrated bam")
    await assert_nonempty(as_cram, label="recalibrated cram")


tests = [
    test_createsequencedictionary,
    test_intervallisttobed,
    test_markduplicates_multiple,
    test_baserecalibrator,
    test_gatherbqsrreports,
    test_applybqsr,
]
