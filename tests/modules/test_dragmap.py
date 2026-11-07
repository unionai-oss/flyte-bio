"""Tests for flyte_bio.modules.dragmap.

The cases and expected values are upstream's dragmap/hashtable and
dragmap/align tests (modules commit 9a1cc9f, DRAGMAP 1.2.1). dragmap/hashtable
snapshots only the output directory's name; here the hash table's main files
are checked instead. dragmap/align snapshots the BAM's name and the log's
lines that start with ``decompHash`` (the hash table loading steps). Every
align case passes the SARS-CoV-2 FASTA as the reference, the human one too, as
upstream; it's unused for BAM output.
"""

import asyncio

from flyte.io import Dir, File

from flyte_bio.modules.dragmap import DragmapAlignResult, dragmap_align, dragmap_hashtable
from tests.framework import env, fixture

SARSCOV2 = "genomics/sarscov2/"
HUMAN = "genomics/homo_sapiens/"
HASH_TABLE_FILES = ["hash_table.cfg", "hash_table.cmp", "reference.bin"]
DECOMPRESSION_LINES = [
    "decompHashTableCtxInit...",
    "decompHashTableHeader...",
    "decompHashTableLiterals...",
    "decompHashTableExtIndex...",
    "decompHashTableAutoHits...",
    "decompHashTableSetFlags...",
]


async def check_alignment(result: DragmapAlignResult, label: str) -> None:
    assert result.alignment.name == "test.bam", f"{label}: alignment is {result.alignment.name}, expected test.bam"
    async with result.log.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    decompression = [line for line in lines if line.startswith("decompHash")]
    assert decompression == DECOMPRESSION_LINES, f"{label}: decompHash log lines {decompression}"


async def reads(genome: str) -> tuple[File, File]:
    return await asyncio.gather(
        fixture(genome + "illumina/fastq/test_1.fastq.gz"), fixture(genome + "illumina/fastq/test_2.fastq.gz")
    )


async def hashtable(genome: str) -> Dir:
    return await dragmap_hashtable(await fixture(genome + "genome/genome.fasta"))


@env.task
async def test_dragmap_hashtable() -> None:
    # upstream case: dragmap/hashtable "sarscov2 - fasta"
    table = await hashtable(SARSCOV2)
    for name in HASH_TABLE_FILES:
        assert await table.get_file(name) is not None, f"dragen-os wrote no {name}"


@env.task
async def test_dragmap_align_sarscov2() -> None:
    # upstream cases: dragmap/align "sarscov2 - fastq, hashtable, fasta, false" / "..., true"
    # and "sarscov2 - [fastq1, fastq2], hashtable, fasta, false" / "..., true"
    table, fasta, (r1, r2) = await asyncio.gather(
        hashtable(SARSCOV2), fixture(SARSCOV2 + "genome/genome.fasta"), reads(SARSCOV2)
    )
    cases = [
        ("single-end, unsorted", r1, None, False),
        ("single-end, sorted", r1, None, True),
        ("paired-end, unsorted", r1, r2, False),
        ("paired-end, sorted", r1, r2, True),
    ]
    results = await asyncio.gather(
        *(dragmap_align(mate_1, mate_2, table, "test", fasta=fasta, sort_bam=sort) for _, mate_1, mate_2, sort in cases)
    )
    for (label, *_), result in zip(cases, results, strict=True):
        await check_alignment(result, label)


@env.task
async def test_dragmap_align_human() -> None:
    # upstream case: dragmap/align "homo_sapiens - [fastq1, fastq2], hashtable, fasta, true"
    table, fasta, (r1, r2) = await asyncio.gather(
        hashtable(HUMAN), fixture(SARSCOV2 + "genome/genome.fasta"), reads(HUMAN)
    )
    result = await dragmap_align(r1, r2, table, "test", fasta=fasta, sort_bam=True)
    await check_alignment(result, "human paired-end, sorted")


tests = [
    test_dragmap_hashtable,
    test_dragmap_align_sarscov2,
    test_dragmap_align_human,
]
