"""Tests for flyte_bio.modules.bwameth.

The cases and expected values are upstream's bwameth/index and bwameth/align
tests (modules commit 9a1cc9f, bwa-meth 0.2.10). As upstream, the alignment
tests build their index from the genome first (with bwa, or bwa-mem2 for the
mem2 case), and compare the BAM by its reads' md5.
"""

from flyte.io import Dir, File

from flyte_bio.modules.bwameth import bwameth_align, bwameth_index
from tests.framework import assert_md5, assert_reads_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"
FASTQ = "genomics/sarscov2/illumina/fastq/"


async def bwameth_check_index(index: Dir, expected: dict[str, str]) -> None:
    files = {f.name: f async for f in index.walk()}
    assert sorted(files) == sorted(expected), f"index files differ: {sorted(files)}"
    for name, md5 in expected.items():
        await assert_md5(files[name], md5, label=name)


async def bwameth_reads(paired: bool) -> list[File]:
    mates = ["1", "2"] if paired else ["1"]
    return [await fixture(f"{FASTQ}test.methylated_{m}.fastq.gz") for m in mates]


C2T = "genome.fasta.bwameth.c2t"
SHARED = {
    C2T: "98039984526a41d04d6bd92fcc040c62",
    f"{C2T}.amb": "249a4195069071ce47cd0bae68abe376",
    f"{C2T}.ann": "46524d4359dcdfb203a235ab3b930dbb",
    f"{C2T}.pac": "4d8e51cb0bbdeaf24576bdf0264d8653",
}


@env.task
async def test_bwameth_index() -> None:
    # upstream case: "sarscov2 genome [fasta]"
    index = await bwameth_index(await fixture(GENOME))
    expected = {
        **SHARED,
        f"{C2T}.bwt": "84f65df7d42dbe84c9ccfaddfdd5ea6b",
        f"{C2T}.sa": "d25f6486f5134f57ed5b258f6fbb8673",
    }
    await bwameth_check_index(index, expected)


@env.task
async def test_bwameth_index_mem2() -> None:
    # upstream case: "sarscov2 genome [fasta] - mem2"
    index = await bwameth_index(await fixture(GENOME), use_mem2=True)
    expected = {
        **SHARED,
        f"{C2T}.0123": "1b92d7081e7b54dafcd97aac45a8076d",
        f"{C2T}.bwt.2bit.64": "52fa820b3a029fe6a9b2e18972f616f2",
    }
    await bwameth_check_index(index, expected)


async def bwameth_aligned(paired: bool, use_mem2: bool = False) -> File:
    fasta = await fixture(GENOME)
    index = await bwameth_index(fasta, use_mem2=use_mem2)
    return await bwameth_align(await bwameth_reads(paired), fasta, index, prefix="test")


@env.task
async def test_bwameth_align_single_end() -> None:
    # upstream case: "sarscov2 methylated single_end [fastq]"
    await assert_reads_md5(await bwameth_aligned(False), "5fdda0ada69daa3956c65a24c9d1c1e0", label="single-end")


@env.task
async def test_bwameth_align_paired_end() -> None:
    # upstream case: "sarscov2 methylated paired_end [fastq]"
    await assert_reads_md5(await bwameth_aligned(True), "cbd5cfa8b446e26987898a45a94770ba", label="paired-end")


@env.task
async def test_bwameth_align_single_end_mem2() -> None:
    # upstream case: "sarscov2 methylated single_end [fastq] - mem2"
    bam = await bwameth_aligned(False, use_mem2=True)
    await assert_reads_md5(bam, "5fdda0ada69daa3956c65a24c9d1c1e0", label="single-end mem2")


tests = [
    test_bwameth_index,
    test_bwameth_index_mem2,
    test_bwameth_align_single_end,
    test_bwameth_align_paired_end,
    test_bwameth_align_single_end_mem2,
]
