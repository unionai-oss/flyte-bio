"""Tests for flyte_bio.modules.art.

The cases and expected values are upstream's art/illumina tests (modules
commit 9a1cc9f, ART 2016.06.05): seeded (``-rs 42``) simulations of HS25
reads from the SARS-CoV-2 genome at 15x, 150 bp. As upstream: the
single-end ``.aln`` md5, and each FASTQ's line count and first record.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.art import art_illumina
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"


async def fastq_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode().splitlines()


async def check_fastq(file: File, count: int, name: str, read: str, quality: str) -> None:
    lines = await fastq_lines(file)
    assert len(lines) == count, f"{file.name}: {len(lines)} lines, expected {count}"
    assert lines[0] == name, f"{file.name}: first read is {lines[0]!r}, expected {name!r}"
    assert lines[1] == read, f"{file.name}: first sequence differs"
    assert lines[3] == quality, f"{file.name}: first quality string differs"


@env.task
async def test_art_illumina_single_end() -> None:
    result = await art_illumina(await fixture(GENOME), "HS25", 15, 150, prefix="test", args="-rs 42")
    assert [f.name for f in result.fastq] == ["test.fq.gz"], [f.name for f in result.fastq]
    assert [f.name for f in result.aln] == ["test.aln"], [f.name for f in result.aln]
    await assert_md5(result.aln[0], "7aa88720911ca1a862e4a24d57435e3c", label="test.aln")
    await check_fastq(
        result.fastq[0],
        11880,
        "@MT192765.1-2970",
        "GAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCTATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCC"
        "ATAATCAAGACTATTCAACCAAGGGTTGAAAAGAAAAAGCTTGATGGC",
        "C1CGGGCGG=GGGJGJG=JJCJJJCJGJJJJGJJJ1JJJGJJGGJGJJGGGJJGJJJCJGGGJ(=G1GGG(CCCGGGCGCGGGCCGCCG8GGGGGGGGGGCC("
        "CGG1GGGG1GGCGCGGCGGGCGCGGGCGGGGCGCCGGGGGCG8GGGG",
    )


@env.task
async def test_art_illumina_paired_end() -> None:
    result = await art_illumina(await fixture(GENOME), "HS25", 15, 150, prefix="test", args="-p -m 250 -s 50 -rs 42")
    assert [f.name for f in result.fastq] == ["test1.fq.gz", "test2.fq.gz"], [f.name for f in result.fastq]
    await check_fastq(
        result.fastq[0],
        5940,
        "@MT192765.1-2970/1",
        "GTTGCGACTACGTGATGAGGAACGAGAAGCGGCTTGACTGCCGCCTCTGCTCCCTTCTGCGTAGAAGCCTTTTGGCAATGTTGTTCCTTGAGGAAGTTGTAGC"
        "ACGATTGCAGCATTGTTAGCAGGATTGCGGGTGCCAATGTGATCTTT",
        "CCCGGGGGGGGGGGCJCJJJJJ1JJJGJJ(JJJGGGJJJJGGJGJJGJJJCCGCJCJC=GCGGG8GJGGGCGG=1G1CJGG=GGCGGCGGGGGCCGG=8GJCC"
        "GCCGCGGGGGGGC=GGG=CCGCGGCGGGCC=GGGGGGCGGCCCGCGG",
    )
    await check_fastq(
        result.fastq[1],
        5940,
        "@MT192765.1-2970/2",
        "GAGCTACCAGACGAATTCGTGGTGGTGACGGTAAAATGAAAGATCTCAGTCCAAGATGGTATTTCTACTACCTAGGAACTGGGCCAGAAGCTGGACTTCCCT"
        "ATGGTGCTAACAAAGACGGCATCATATGGGTTGCAACTGAGGGAGCCT",
        "1CCGGGGGGGCGGGJJJJJGJG1GJGJ=JGJJCJJJCJJJGJ1GJGJGGJJC8JGG=JGGCGCCGCGGCCCGJCGGGGCJC==CGG88CGGGCGGGCCCGGC8"
        "JJJJCGGCGGCCGGGG=CGGGGCCGGCCCGCGCCCCG8GGGG=GGGC",
    )


tests = [test_art_illumina_single_end, test_art_illumina_paired_end]
