"""Tests for flyte_bio.modules.bowtie (Bowtie version 1).

The cases and expected values are upstream's bowtie/build and bowtie/align
tests (modules commit 9a1cc9f, Bowtie 1.3.1). bowtie/build snapshots the index
files' md5s. bowtie/align snapshots, for single-end reads, the alignment's read
md5 (nft-bam ``getReadsMD5``), the decompressed unaligned reads' md5 and the
log's md5; for paired-end reads, which upstream calls unstable, the first two
records' sequences (nft-bam ``getReads(2)``) and the log's md5. Both align
against the index built with prefix ``sarscov2``.
"""

import asyncio

from flyte.io import Dir, File

from flyte_bio.modules.bowtie import bowtie1_align, bowtie1_build
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_gunzipped_md5, assert_md5, assert_reads_md5, env, fixture

FASTA = "genomics/sarscov2/genome/genome.fasta"
FASTQ = "genomics/sarscov2/illumina/fastq/"
PAIRED_END_FIRST_READS = [
    (
        "ATGTGTACATTGGCGACCCTGCTCAATTACCTGCACCACGCACATTGCTAACTAAGGGCACACTAGAACCAGAATATTTCAATTCAGTGTGTAGACTTATGAAAACTA"
        "TAGGTCCAGACATGTTCCTCGGAACTTGTCGGCGTTGTCCTG"
    ),
    (
        "ACGCACATTGCTAACTAAGGGCACACTAGAACCAGAATATTTCAATTCAGTGTGTAGACTTATGAAAACTATAGGTCCAGACATGTTCCTCGGAACTTGTCGGCGTTG"
        "TCCTGCTGAAATTGTTGACACTGTGAGTGCTTTGGTTTATGA"
    ),
]


async def first_reads(bam: File, count: int) -> list[str]:
    """nft-bam's ``getReads(count)``: the SEQ of the first ``count`` records, in file order."""
    sam = await samtools_view(bam=bam)
    async with sam.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    return [line.split("\t")[9] for line in lines[:count]]


async def sarscov2_index() -> Dir:
    return await bowtie1_build(await fixture(FASTA), "sarscov2")


@env.task
async def test_bowtie1_build() -> None:
    # upstream case: bowtie/build "sarscov2 - fasta"
    index = await sarscov2_index()
    expected = {
        "sarscov2.1.ebwt": "d9b76ecf9fd0413240173273b38d8199",
        "sarscov2.2.ebwt": "02b44af9f94c62ecd3c583048e25d4cf",
        "sarscov2.3.ebwt": "4ed93abba181d8dfab2e303e33114777",
        "sarscov2.4.ebwt": "c25be5f8b0378abf7a58c8a880b87626",
        "sarscov2.rev.1.ebwt": "b37aaf11853e65a3b13561f27a912b06",
        "sarscov2.rev.2.ebwt": "9e6b0c4c1ddb99ae71ff8a4fe5ec6459",
    }
    for name, md5 in expected.items():
        file = await index.get_file(name)
        assert file is not None, f"bowtie-build wrote no {name}"
        await assert_md5(file, md5, label=name)


@env.task
async def test_bowtie1_align_single_end() -> None:
    # upstream case: bowtie/align "sarscov2 - single_end" (save_unaligned true)
    index, reads = await asyncio.gather(sarscov2_index(), fixture(FASTQ + "test_1.fastq.gz"))
    result = await bowtie1_align(reads, None, index, "test", save_unaligned=True)
    assert result.bam.name == "test.bam", f"alignment is {result.bam.name}, expected test.bam"
    await assert_reads_md5(result.bam, "7bdcfc6f54ae6e8f4570395cc85db9a3", label="single-end")
    await assert_md5(result.log, "4b9140ceadb8a18ae9330885370f8a0b", label="single-end test.out")
    assert [f.name for f in result.unaligned] == ["test.unmapped.fastq.gz"], f"unaligned {result.unaligned}"
    await assert_gunzipped_md5(result.unaligned[0], "5729a694abd09657da3b9101861090c4", label="unaligned")


@env.task
async def test_bowtie1_align_paired_end() -> None:
    # upstream case: bowtie/align "sarscov2 - paired_end" (save_unaligned false)
    index, r1, r2 = await asyncio.gather(
        sarscov2_index(), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    result = await bowtie1_align(r1, r2, index, "test")
    assert result.bam.name == "test.bam", f"alignment is {result.bam.name}, expected test.bam"
    reads = await first_reads(result.bam, 2)
    assert reads == PAIRED_END_FIRST_READS, f"paired-end first reads {reads}"
    await assert_md5(result.log, "5e13272d112cef8faeedcdbd7c602de0", label="paired-end test.out")


tests = [
    test_bowtie1_build,
    test_bowtie1_align_single_end,
    test_bowtie1_align_paired_end,
]
