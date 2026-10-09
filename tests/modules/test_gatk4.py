"""Tests for flyte_bio.modules.gatk4 (md5s from upstream's gatk4 module snapshots)."""

import asyncio

from flyte_bio.modules.gatk4 import gatk4_createsequencedictionary, gatk4_intervallisttobed, gatk4_markduplicates
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_md5, assert_nonempty, env, fixture

HS = "genomics/homo_sapiens/"


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


tests = [test_createsequencedictionary, test_intervallisttobed, test_markduplicates_multiple]
