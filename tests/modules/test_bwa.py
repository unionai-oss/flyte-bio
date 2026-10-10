"""Tests for flyte_bio.modules.bwa (md5s from upstream's bwa/index and bwa/mem snapshots).

bwa/mem's snapshot hashes the alignments' read sequences (nft-bam
``getReadsMD5``), in file order.
"""

import asyncio

from flyte_bio.modules.bwa import bwa_index, bwa_mem
from tests.framework import assert_md5, assert_reads_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_bwa_index() -> None:
    # upstream case: bwa/index "BWA index"
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    index = await bwa_index(fasta)
    expected = {
        "genome.amb": "3a68b8b2287e07dd3f5f95f4344ba76e",
        "genome.ann": "c32e11f6c859f166c7525a9c1d583567",
        "genome.bwt": "0469c30a1e239dd08f68afe66fde99da",
        "genome.pac": "983e3d2cd6f36e2546e6d25a0da78d66",
        "genome.sa": "ab3952cabf026b48cd3eb5bccbb636d1",
    }
    for name, md5 in expected.items():
        f = await index.get_file(name)
        assert f is not None, f"bwa index wrote no {name}"
        await assert_md5(f, md5, label=name)


@env.task
async def test_bwa_mem() -> None:
    # upstream cases: bwa/mem "Single-End", "Single-End Sort", "Paired-End", "Paired-End Sort"
    fasta, r1, r2 = await asyncio.gather(
        fixture("genomics/sarscov2/genome/genome.fasta"),
        fixture(FASTQ + "test_1.fastq.gz"),
        fixture(FASTQ + "test_2.fastq.gz"),
    )
    index = await bwa_index(fasta)
    cases = [
        ("single-end", r1, None, False, "798439cbd7fd81cbcc5078022dc5479d"),
        ("single-end sorted", r1, None, True, "94fcf617f5b994584c4e8d4044e16b4f"),
        ("paired-end", r1, r2, False, "57aeef88ed701a8ebc8e2f0a381b2a6"),
        ("paired-end sorted", r1, r2, True, "af8628d9df18b2d3d4f6fd47ef2bb872"),
    ]
    bams = await asyncio.gather(*(bwa_mem(a, b, index, prefix="test", sort=sort) for _, a, b, sort, _ in cases))
    for (label, *_, md5), bam in zip(cases, bams):
        assert bam.path.endswith("/test.bam"), bam.path
        await assert_reads_md5(bam, md5, label=label)


tests = [test_bwa_index, test_bwa_mem]
