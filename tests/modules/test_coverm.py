"""Tests for flyte_bio.modules.coverm.

The cases and expected md5s are upstream's coverm tests (modules commit
9a1cc9f, CoverM 0.7.0). As upstream, the BAMs CoverM writes are checked by
their reads' md5 (``getReadsMD5``). The module's 4 CPUs match upstream's test
CPUs. As upstream's config for the directory case, it reads the genomes with
``--genome-fasta-extension fasta``; upstream's "invalid ref_mode" case, which
expects the module to fail, checks that the wrapper raises.
"""

from flyte.io import File

from flyte_bio.modules.coverm import CovermResult, coverm_contig, coverm_genome
from tests.framework import assert_md5, assert_reads_md5, env, fixture, fixture_dir

FRAGILIS = "genomics/prokaryotes/bacteroides_fragilis/"
BAM_READS = ["edb6d55735486ce7c8821fcb0e65610d", "910bd83216afed9a16f12fdfc074366b"]


async def fragilis_reads() -> list[File]:
    names = ["test1_1", "test1_2", "test2_1", "test2_2"]
    return [await fixture(f"{FRAGILIS}illumina/fastq/{n}.fastq.gz") for n in names]


async def fragilis_bams() -> list[File]:
    return [await fixture(f"{FRAGILIS}illumina/bam/{n}.sorted.bam") for n in ["test1", "test2"]]


async def fragilis_genome() -> File:
    return await fixture(FRAGILIS + "genome/genome.fna.gz")


async def check(result: CovermResult, coverage: str, bam_reads: list[str]) -> None:
    await assert_md5(result.coverage, coverage, label="depth")
    assert len(result.bams) == len(bam_reads), f"CoverM wrote {len(result.bams)} BAMs, expected {len(bam_reads)}"
    for bam, md5 in zip(result.bams, bam_reads, strict=True):
        await assert_reads_md5(bam, md5, label=bam.name)


@env.task
async def test_coverm_contig_fastq() -> None:
    result = await coverm_contig(await fragilis_reads(), await fragilis_genome(), enable_bam_output=True, prefix="test")
    await check(result, "2d68cbaf0dc3c88a024197386a68d71d", BAM_READS)


@env.task
async def test_coverm_contig_bam() -> None:
    result = await coverm_contig(await fragilis_bams(), bam_input=True, prefix="test")
    await check(result, "8e97198137f9ab7bdca3124f3e544438", [])


@env.task
async def test_coverm_genome_fastq() -> None:
    result = await coverm_genome(
        await fragilis_reads(), await fragilis_genome(), ref_mode="file", enable_bam_output=True, prefix="test"
    )
    await check(result, "c16918076fbb85545cd0f29a228cb54f", BAM_READS)


@env.task
async def test_coverm_genome_bam() -> None:
    result = await coverm_genome(
        await fragilis_bams(), await fragilis_genome(), bam_input=True, ref_mode="", prefix="test"
    )
    await check(result, "60ac75e5dce5da12263fbeea80f9c7fa", [])


@env.task
async def test_coverm_genome_directory() -> None:
    genomes = await fixture_dir("genomics/homo_sapiens/genome/chr21/sequence/chromosomes.tar.gz", "chromosomes")
    reads = [await fixture(f"genomics/homo_sapiens/illumina/fastq/test_{n}.fastq.gz") for n in (1, 2)]
    result = await coverm_genome(reads, genomes, ref_mode="dir", prefix="test", args="--genome-fasta-extension fasta")
    await check(result, "70b9ad4a4d90d9c6eb0ef4b52adaf636", [])


@env.task
async def test_coverm_genome_invalid_ref_mode() -> None:
    reads = (await fragilis_reads())[:2]
    try:
        await coverm_genome(reads, await fragilis_genome(), ref_mode="not_a_valid_mode", prefix="test")
    except ValueError as error:
        assert "Invalid ref_mode" in str(error), f"unexpected error {error}"
    else:
        raise AssertionError("coverm_genome accepted ref_mode 'not_a_valid_mode'")


tests = [
    test_coverm_contig_fastq,
    test_coverm_contig_bam,
    test_coverm_genome_fastq,
    test_coverm_genome_bam,
    test_coverm_genome_directory,
    test_coverm_genome_invalid_ref_mode,
]
