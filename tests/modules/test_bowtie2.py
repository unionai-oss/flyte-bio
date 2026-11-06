"""Tests for flyte_bio.modules.bowtie2.

The cases and expected values are upstream's bowtie2/build and bowtie2/align
tests (modules commit 9a1cc9f, bowtie2 2.5.5). bowtie2/build snapshots the
index files' md5s. bowtie2/align snapshots the alignment's name, nft-bam's
``getStatistics`` (read lengths, mapping qualities, read and duplicate counts,
recomputed here by :func:`alignment_stats`), the log's md5 and the decompressed
unaligned reads' md5s. Its tests build the index with ``--seed 1`` and align
with ``--seed 1``, except the CRAM case, which sets neither.
"""

import asyncio

from flyte.io import Dir, File

from flyte_bio.modules.bowtie2 import Bowtie2AlignResult, bowtie2_align, bowtie2_build
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTA = "genomics/sarscov2/genome/genome.fasta"
FASTQ = "genomics/sarscov2/illumina/fastq/"
SINGLE_END_LOG = "2318a3c202cbd5893b9a975c162f4fdd"
PAIRED_END_LOG = "169b09bba82894cb730cc16512541820"


def stats(min_read_length: int, read_count: int) -> dict[str, int]:
    """Upstream's snapshot statistics; only the shortest read and the read count differ between cases."""
    return {
        "maxReadLength": 151,
        "minReadLength": min_read_length,
        "meanReadLength": 138,
        "maxQuality": 42,
        "minQuality": 0,
        "meanQuality": 41,
        "readCount": read_count,
        "duplicateReadCount": 0,
    }


SINGLE_END_STATS = stats(72, 100)
PAIRED_END_STATS = stats(53, 200)


async def alignment_stats(alignment: File, fasta: File | None = None) -> dict[str, int]:
    """nft-bam's ``getStatistics``: over every record, SEQ lengths and MAPQs (means truncated), counts."""
    sam = await samtools_view(bam=alignment, fasta=fasta)
    async with sam.open("rb") as fh:
        records = [line.split("\t") for line in bytes(await fh.read()).decode().splitlines()]
    lengths = [0 if r[9] == "*" else len(r[9]) for r in records]
    qualities = [int(r[4]) for r in records]
    return {
        "maxReadLength": max(lengths),
        "minReadLength": min(lengths),
        "meanReadLength": sum(lengths) // len(records),
        "maxQuality": max(qualities),
        "minQuality": min(qualities),
        "meanQuality": sum(qualities) // len(records),
        "readCount": len(records),
        "duplicateReadCount": sum(1 for r in records if int(r[1]) & 0x400),
    }


async def check_alignment(
    result: Bowtie2AlignResult,
    name: str,
    expected: dict[str, int],
    label: str,
    log_md5: str | None = None,
    unaligned_md5s: tuple[str, ...] = (),
    fasta: File | None = None,
) -> None:
    assert result.alignment.name == name, f"{label}: alignment is {result.alignment.name}, expected {name}"
    actual = await alignment_stats(result.alignment, fasta)
    assert actual == expected, f"{label}: statistics {actual}, expected {expected}"
    if log_md5 is not None:
        await assert_md5(result.log, log_md5, label=f"{label} log")
    assert len(result.unaligned) == len(unaligned_md5s), f"{label}: {len(result.unaligned)} unaligned files"
    for file, md5 in zip(result.unaligned, unaligned_md5s, strict=True):
        await assert_gunzipped_md5(file, md5, label=f"{label} {file.name}")


async def seeded_index(args: str = "--seed 1") -> Dir:
    return await bowtie2_build(await fixture(FASTA), args=args)


@env.task
async def test_bowtie2_build() -> None:
    # upstream case: bowtie2/build "Should run without failures"
    index = await bowtie2_build(await fixture(FASTA))
    expected = {
        "genome.1.bt2": "cbe3d0bbea55bc57c99b4bfa25b5fbdf",
        "genome.2.bt2": "47b153cd1319abc88dda532462651fcf",
        "genome.3.bt2": "4ed93abba181d8dfab2e303e33114777",
        "genome.4.bt2": "c25be5f8b0378abf7a58c8a880b87626",
        "genome.rev.1.bt2": "52be6950579598a990570fbcf5372184",
        "genome.rev.2.bt2": "e3b4ef343dea4dd571642010a7d09597",
    }
    for name, md5 in expected.items():
        file = await index.get_file(name)
        assert file is not None, f"bowtie2-build wrote no {name}"
        await assert_md5(file, md5, label=name)


@env.task
async def test_bowtie2_align_single_end() -> None:
    # upstream cases: bowtie2/align "fastq, index, fasta_fai, false, true - bam", "... - sam", "... - sam2"
    # and "fastq, index, fasta_fai, true, true - bam"
    index, fasta, reads = await asyncio.gather(seeded_index(), fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"))
    cases = [
        ("bam", "test.bam", False, ""),
        ("sam", "test.sam", False, "--output-fmt SAM"),
        ("sam2", "test.sam", False, "-O SAM"),
        ("bam, unaligned", "test.bam", True, ""),
    ]
    results = await asyncio.gather(
        *(
            bowtie2_align(reads, None, index, "test", fasta=fasta, save_unaligned=save, args="--seed 1", args2=args2)
            for _, _, save, args2 in cases
        )
    )
    for (label, name, save, _), result in zip(cases, results, strict=True):
        unaligned = ("d17d9a4ff728b6a6afa6203a746b7a01",) if save else ()
        await check_alignment(result, name, SINGLE_END_STATS, f"single-end {label}", SINGLE_END_LOG, unaligned)


@env.task
async def test_bowtie2_align_paired_end() -> None:
    # upstream cases: bowtie2/align "[fastq1, fastq2], index, fasta_fai, false, true - bam"
    # and "[fastq1, fastq2], index, fasta_fai, true, true - bam"
    index, fasta, r1, r2 = await asyncio.gather(
        seeded_index(), fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    plain, unaligned = await asyncio.gather(
        bowtie2_align(r1, r2, index, "test", fasta=fasta, args="--seed 1"),
        bowtie2_align(r1, r2, index, "test", fasta=fasta, save_unaligned=True, args="--seed 1"),
    )
    await check_alignment(plain, "test.bam", PAIRED_END_STATS, "paired-end", PAIRED_END_LOG)
    await check_alignment(
        unaligned,
        "test.bam",
        PAIRED_END_STATS,
        "paired-end unaligned",
        PAIRED_END_LOG,
        ("84b4ab15e6ffc9f29d0a3f1b79003744", "40bbe14edd48a98cca6f53574b2198fd"),
    )


@env.task
async def test_bowtie2_align_large_index() -> None:
    # upstream cases: bowtie2/align "fastq, large_index, fasta_fai, false, true - bam"
    # and "[fastq1, fastq2], large_index, fasta_fai, false, true - bam"
    index, fasta, r1, r2 = await asyncio.gather(
        seeded_index("--large-index --seed 1"),
        fixture(FASTA),
        fixture(FASTQ + "test_1.fastq.gz"),
        fixture(FASTQ + "test_2.fastq.gz"),
    )
    assert await index.get_file("genome.rev.1.bt2l") is not None, "--large-index wrote no genome.rev.1.bt2l"
    single, paired = await asyncio.gather(
        bowtie2_align(r1, None, index, "test", fasta=fasta, args="--seed 1"),
        bowtie2_align(r1, r2, index, "test", fasta=fasta, args="--seed 1"),
    )
    await check_alignment(single, "test.bam", SINGLE_END_STATS, "large index single-end", SINGLE_END_LOG)
    await check_alignment(paired, "test.bam", PAIRED_END_STATS, "large index paired-end", PAIRED_END_LOG)


@env.task
async def test_bowtie2_align_cram() -> None:
    # upstream case: bowtie2/align "[fastq1, fastq2], index, fasta_fai, false, true - cram"
    fasta, r1, r2 = await asyncio.gather(
        fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    index = await bowtie2_build(fasta)
    result = await bowtie2_align(r1, r2, index, "test", fasta=fasta, args2="--output-fmt cram --write-index")
    assert result.alignment_index is not None and result.alignment_index.name == "test.cram.crai", (
        f"cram index is {result.alignment_index and result.alignment_index.name}, expected test.cram.crai"
    )
    await check_alignment(result, "test.cram", PAIRED_END_STATS, "paired-end cram", fasta=fasta)


tests = [
    test_bowtie2_build,
    test_bowtie2_align_single_end,
    test_bowtie2_align_paired_end,
    test_bowtie2_align_large_index,
    test_bowtie2_align_cram,
]
