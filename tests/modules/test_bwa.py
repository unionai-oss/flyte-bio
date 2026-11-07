"""Tests for flyte_bio.modules.bwa (md5s from upstream's bwa/index, bwa/mem, bwa/samse and bwa/sampe snapshots).

bwa/mem's, bwa/samse's and bwa/sampe's snapshots hash the alignments' read
sequences (nft-bam ``getReadsMD5``), in file order. bwa/aln's snapshot holds
only the ``.sai`` names (modules commit 9a1cc9f, bwa 0.7.19).
"""

import asyncio

from flyte.io import Dir, File

from flyte_bio.modules.bwa import bwa_aln, bwa_index, bwa_mem, bwa_sampe, bwa_samse
from tests.framework import assert_md5, assert_nonempty, assert_reads_md5, env, fixture

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


async def sarscov2_reads() -> tuple[File, File, Dir]:
    fasta, r1, r2 = await asyncio.gather(
        fixture("genomics/sarscov2/genome/genome.fasta"),
        fixture(FASTQ + "test_1.fastq.gz"),
        fixture(FASTQ + "test_2.fastq.gz"),
    )
    return r1, r2, await bwa_index(fasta)


@env.task
async def test_bwa_aln_single_end() -> None:
    # upstream case: bwa/aln "Single-End"
    r1, _, index = await sarscov2_reads()
    sais = await bwa_aln(r1, None, index, prefix="test")
    assert [s.name for s in sais] == ["test.sai"], [s.name for s in sais]
    await assert_nonempty(sais[0], label="test.sai")


@env.task
async def test_bwa_aln_paired_end() -> None:
    # upstream case: bwa/aln "Paired-End"
    r1, r2, index = await sarscov2_reads()
    sais = await bwa_aln(r1, r2, index, prefix="test")
    assert [s.name for s in sais] == ["test.1.sai", "test.2.sai"], [s.name for s in sais]
    for sai in sais:
        await assert_nonempty(sai, label=sai.name)


@env.task
async def test_bwa_samse() -> None:
    # upstream case: bwa/samse "SamSE"
    r1, _, index = await sarscov2_reads()
    (sai,) = await bwa_aln(r1, None, index, prefix="test")
    bam = await bwa_samse(r1, sai, index, prefix="test")
    assert bam.name == "test.bam", bam.name
    await assert_reads_md5(bam, "1ca0f48b86ee8cf2417f72d6bf95939a", label="test.bam")


@env.task
async def test_bwa_sampe() -> None:
    # upstream case: bwa/sampe "SamPE"
    r1, r2, index = await sarscov2_reads()
    sais = await bwa_aln(r1, r2, index, prefix="test")
    bam = await bwa_sampe(r1, r2, sais, index, prefix="test")
    assert bam.name == "test.bam", bam.name
    await assert_reads_md5(bam, "ae4005223c9eb0f2b11bcd984e088f9c", label="test.bam")


tests = [
    test_bwa_index,
    test_bwa_mem,
    test_bwa_aln_single_end,
    test_bwa_aln_paired_end,
    test_bwa_samse,
    test_bwa_sampe,
]
