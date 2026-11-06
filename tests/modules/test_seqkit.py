"""Tests for flyte_bio.modules.seqkit.

The cases and expected md5s are upstream's seqkit/stats and seqkit/replace
tests (modules commit 9a1cc9f, seqkit 2.13.0), with each case's ``ext.args``.
Upstream's replace patterns are shell-quoted (``-p 'A'``); here they are bare.
nf-test hashes ``.gz`` outputs decompressed.
"""

from flyte.io import File

from flyte_bio.modules.seqkit import seqkit_replace, seqkit_stats
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
GENOME = "genomics/sarscov2/genome/"


def name_of(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


@env.task
async def test_seqkit_stats_single_end() -> None:
    stats = await seqkit_stats([await fixture(FASTQ + "test_1.fastq.gz")], prefix="test")
    await assert_md5(stats, "76d5ab444de880ef443f358841e44cf9", label="test.tsv")


@env.task
async def test_seqkit_stats_paired_end() -> None:
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    stats = await seqkit_stats(reads, prefix="test")
    await assert_md5(stats, "2d96d9f7d85456090e138b81bf97bcd5", label="test.tsv")


@env.task
async def test_seqkit_stats_nanopore() -> None:
    stats = await seqkit_stats([await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz")], prefix="test")
    await assert_md5(stats, "388047b0263b8d0a3041d6d80b316782", label="test.tsv")


@env.task
async def test_seqkit_stats_genome_fasta() -> None:
    stats = await seqkit_stats([await fixture(GENOME + "genome.fasta")], prefix="test")
    await assert_md5(stats, "6f6f8b42abd183c6ab90346a311ca522", label="test.tsv")


@env.task
async def test_seqkit_stats_transcriptome_fasta() -> None:
    stats = await seqkit_stats([await fixture(GENOME + "transcriptome.fasta")], prefix="test")
    await assert_md5(stats, "c79d860b2b0379810e4df3e8359fc8ed", label="test.tsv")


@env.task
async def test_seqkit_replace_fasta_gz() -> None:
    out = await seqkit_replace(await fixture(GENOME + "proteome.fasta.gz"), prefix="test", args="-s -p A -r N")
    assert name_of(out) == "test.fasta.gz", f"expected test.fasta.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, "c40eaff961f6f2a48bb7e8fd156ed5d7", label="test.fasta.gz")


@env.task
async def test_seqkit_replace_uncompressed() -> None:
    out = await seqkit_replace(await fixture(GENOME + "genome.fasta"), prefix="test", args="-s -p A -r N")
    assert name_of(out) == "test.fasta", f"expected test.fasta, got {name_of(out)}"
    await assert_md5(out, "b1518908253a4997fcad98270751112e", label="test.fasta")


@env.task
async def test_seqkit_replace_custom_extension() -> None:
    genome = await fixture(GENOME + "genome.fasta")
    out = await seqkit_replace(genome, prefix="test", out_ext="fasta", args="-s -p T -r N")
    assert name_of(out) == "test.fasta", f"expected test.fasta, got {name_of(out)}"
    await assert_md5(out, "05d3294a62c72f5489f067c1da3c2f6c", label="test.fasta")


tests = [
    test_seqkit_stats_single_end,
    test_seqkit_stats_paired_end,
    test_seqkit_stats_nanopore,
    test_seqkit_stats_genome_fasta,
    test_seqkit_stats_transcriptome_fasta,
    test_seqkit_replace_fasta_gz,
    test_seqkit_replace_uncompressed,
    test_seqkit_replace_custom_extension,
]
