"""Tests for flyte_bio.modules.hisat2.

The cases and expected values are upstream's hisat2/extractsplicesites,
hisat2/build and hisat2/align tests (modules commit 9a1cc9f, HISAT2 2.2.3).
Upstream marks the splice sites, the index and the BAMs unstable, so their
snapshots hold only the file names; hisat2/align also snapshots the summary
log's md5 and that no unaligned reads were written. Upstream builds the index
with a memory threshold of 3 GB (splice-aware) for the cases with splice
sites, and with ``false`` (never splice-aware) otherwise.
"""

import asyncio

from flyte.io import Dir, File

from flyte_bio.modules.hisat2 import Hisat2AlignResult, hisat2_align, hisat2_build, hisat2_extractsplicesites
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/"
FASTQ = "genomics/sarscov2/illumina/fastq/"
INDEX_FILES = [f"genome.{n}.ht2" for n in range(1, 9)]
SINGLE_END_SUMMARY = "7b8a9e61b7646da1089b041333c41a87"
PAIRED_END_SUMMARY = "9839b31db795958cc4b70711a3414e9c"


async def check_index(index: Dir, label: str) -> None:
    names = sorted(f.name for f in await index.list_files())
    assert names == INDEX_FILES, f"{label}: index files {names}, expected {INDEX_FILES}"


async def check_alignment(result: Hisat2AlignResult, summary_md5: str, label: str) -> None:
    assert result.bam.name == "test.bam", f"{label}: alignment is {result.bam.name}, expected test.bam"
    assert result.summary.name == "test.hisat2.summary.log", f"{label}: summary is {result.summary.name}"
    await assert_md5(result.summary, summary_md5, label=f"{label} summary")
    names = [f.name for f in await result.results.list_files()]
    unaligned = [name for name in names if name.endswith("fastq.gz")]
    assert not unaligned, f"{label}: unexpected unaligned reads {unaligned}"


async def splice_aware_index() -> tuple[Dir, File]:
    """Upstream's setup for the cases with splice sites: extract them, then build with a 3 GB threshold."""
    fasta, gtf = await asyncio.gather(fixture(GENOME + "genome.fasta"), fixture(GENOME + "genome.gtf"))
    splicesites = await hisat2_extractsplicesites(gtf)
    index = await hisat2_build(fasta, gtf=gtf, splicesites=splicesites, splice_aware_memory_gb=3)
    return index, splicesites


@env.task
async def test_hisat2_extractsplicesites() -> None:
    # upstream case: hisat2/extractsplicesites "Should run without failures"
    splicesites = await hisat2_extractsplicesites(await fixture(GENOME + "genome.gtf"))
    assert splicesites.name == "genome.splice_sites.txt", f"splice sites are {splicesites.name}"


@env.task
async def test_hisat2_build() -> None:
    # upstream cases: hisat2/build "Should run without failures" (threshold false)
    # and "Should use splice sites when memory threshold is met" (threshold 3 GB)
    fasta, gtf = await asyncio.gather(fixture(GENOME + "genome.fasta"), fixture(GENOME + "genome.gtf"))
    splicesites = await hisat2_extractsplicesites(gtf)
    plain, splice_aware = await asyncio.gather(
        hisat2_build(fasta, gtf=gtf, splicesites=splicesites),
        hisat2_build(fasta, gtf=gtf, splicesites=splicesites, splice_aware_memory_gb=3),
    )
    await check_index(plain, "no threshold")
    await check_index(splice_aware, "3 GB threshold")


@env.task
async def test_hisat2_align_splice_sites() -> None:
    # upstream cases: hisat2/align "Single-End" and "Paired-End"
    (index, splicesites), r1, r2 = await asyncio.gather(
        splice_aware_index(), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    single, paired = await asyncio.gather(
        hisat2_align(r1, None, index, "test", splicesites=splicesites),
        hisat2_align(r1, r2, index, "test", splicesites=splicesites),
    )
    await check_alignment(single, SINGLE_END_SUMMARY, "single-end")
    await check_alignment(paired, PAIRED_END_SUMMARY, "paired-end")


@env.task
async def test_hisat2_align_no_splice_sites() -> None:
    # upstream cases: hisat2/align "Single-End No Splice Sites" and "Paired-End No Splice Sites"
    fasta, r1, r2 = await asyncio.gather(
        fixture(GENOME + "genome.fasta"), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    index = await hisat2_build(fasta)
    single, paired = await asyncio.gather(
        hisat2_align(r1, None, index, "test"),
        hisat2_align(r1, r2, index, "test"),
    )
    await check_alignment(single, SINGLE_END_SUMMARY, "single-end, no splice sites")
    await check_alignment(paired, PAIRED_END_SUMMARY, "paired-end, no splice sites")


tests = [
    test_hisat2_extractsplicesites,
    test_hisat2_build,
    test_hisat2_align_splice_sites,
    test_hisat2_align_no_splice_sites,
]
