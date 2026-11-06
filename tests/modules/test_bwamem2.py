"""Tests for flyte_bio.modules.bwamem2 (md5s from upstream's bwamem2/index and bwamem2/mem snapshots).

bwamem2/mem's snapshot hashes the alignments' read sequences (nft-bam
``getReadsMD5``), in file order, and the BAM header (``getHeaderMD5``). The
header records the command line (paths, thread count), so it is checked by
content instead of md5.
"""

import asyncio

from flyte.io import File

from flyte_bio.modules.bwamem2 import bwamem2_index, bwamem2_mem
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_md5, assert_reads_md5, env, fixture

GENOME = "genomics/sarscov2/genome/"
FASTQ = "genomics/sarscov2/illumina/fastq/"

# The same index files, whether built from genome.fasta or genome.fasta.gz.
INDEX_MD5S = {
    "0123": "b02870de80106104abcb03cd9463e7d8",
    "amb": "3a68b8b2287e07dd3f5f95f4344ba76e",
    "ann": "c32e11f6c859f166c7525a9c1d583567",
    "bwt.2bit.64": "d097a1b82dee375d41a1ea69895a9216",
    "pac": "983e3d2cd6f36e2546e6d25a0da78d66",
}


@env.task
async def test_bwamem2_index() -> None:
    # upstream cases: bwamem2/index "fasta" and "fasta.gz"
    names = ["genome.fasta", "genome.fasta.gz"]
    fastas = await asyncio.gather(*(fixture(GENOME + name) for name in names))
    indexes = await asyncio.gather(*(bwamem2_index(fasta) for fasta in fastas))
    for name, index in zip(names, indexes):
        for suffix, md5 in INDEX_MD5S.items():
            f = await index.get_file(f"{name}.{suffix}")
            assert f is not None, f"bwa-mem2 index wrote no {name}.{suffix}"
            await assert_md5(f, md5, label=f"{name}.{suffix}")


async def header_lines(bam: File) -> list[str]:
    header = await samtools_view(bam=bam, args="-H")
    async with header.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


@env.task
async def test_bwamem2_mem() -> None:
    # upstream cases: bwamem2/mem "sarscov2 - fastq, index, fasta, false" / "..., true",
    # "sarscov2 - [fastq1, fastq2], index, fasta, false" / "..., true"
    fasta, r1, r2 = await asyncio.gather(
        fixture(GENOME + "genome.fasta"),
        fixture(FASTQ + "test_1.fastq.gz"),
        fixture(FASTQ + "test_2.fastq.gz"),
    )
    index = await bwamem2_index(fasta)
    cases = [
        ("single-end", r1, None, False, "798439cbd7fd81cbcc5078022dc5479d"),
        ("single-end sorted", r1, None, True, "94fcf617f5b994584c4e8d4044e16b4f"),
        ("paired-end", r1, r2, False, "57aeef88ed701a8ebc8e2f0a381b2a6"),
        ("paired-end sorted", r1, r2, True, "af8628d9df18b2d3d4f6fd47ef2bb872"),
    ]
    bams = await asyncio.gather(
        *(bwamem2_mem(a, b, index, prefix="test", fasta=fasta, sort=sort) for _, a, b, sort, _ in cases)
    )
    for (label, _, _, sort, md5), bam in zip(cases, bams):
        assert bam.path.endswith("/test.bam"), bam.path
        await assert_reads_md5(bam, md5, label=label)
        header = await header_lines(bam)
        assert "@SQ\tSN:MT192765.1\tLN:29829" in header, f"{label}: {header}"
        assert any(line.startswith("@PG\tID:bwa-mem2") for line in header), f"{label}: {header}"
        sorted_header = any(line.startswith("@HD") and "SO:coordinate" in line for line in header)
        assert sorted_header == sort, f"{label}: {header}"


tests = [test_bwamem2_index, test_bwamem2_mem]
