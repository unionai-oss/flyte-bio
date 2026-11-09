"""Tests for flyte_bio.modules.ivar.

The cases and expected md5s are upstream's ivar/consensus and ivar/variants
tests (modules commit 9a1cc9f, iVar 1.4.4). The consensus tests use the test
directory's samtools mpileup options (``-aa -A -d 0 -Q 0``); the variants
tests use none. Upstream's consensus-with-pileup case checks only the names
of the qualities and the pileup; so does its port here.

Upstream's consensus tests pile up ``test.paired_end.bam``, which isn't
sorted: samtools mpileup stops at its first out-of-order read, and the
snapshot is of iVar's consensus of that empty pileup (a header and blank
sequence). The module, as upstream, takes iVar's status for the pipeline's,
so these tests reproduce it.
"""

from flyte_bio.modules.ivar import ivar_consensus, ivar_variants
from tests.framework import assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/test.paired_end.bam"
SORTED_BAM = "genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"
GENOME = "genomics/sarscov2/genome/genome.fasta"
CONSENSUS_MPILEUP_ARGS = "-aa -A -d 0 -Q 0"
CONSENSUS_MD5 = "9e21a64818f4302b4dece5480fa5e8b8"
VARIANTS_MD5 = "76c93420121bb891fa275f716142a157"
VARIANTS_GFF_MD5 = "2cf6c48b55e0e81155ff310d48de94df"
MPILEUP_MD5 = "958e6bead4103d72026f80153b6b5150"


@env.task
async def test_ivar_consensus() -> None:
    # upstream case: ivar/consensus "test-ivar-consensus"
    result = await ivar_consensus(
        await fixture(BAM), await fixture(GENOME), prefix="test", mpileup_args=CONSENSUS_MPILEUP_ARGS
    )
    await assert_md5(result.fasta, CONSENSUS_MD5, label="test.fa")
    await assert_md5(result.qual, "68b329da9893e34099c7d8ad5cb9c940", label="test.qual.txt")
    assert result.mpileup is None


@env.task
async def test_ivar_consensus_mpileup() -> None:
    # upstream case: ivar/consensus "test-ivar-consensus-mpileup"
    result = await ivar_consensus(
        await fixture(BAM),
        await fixture(GENOME),
        save_mpileup=True,
        prefix="test",
        mpileup_args=CONSENSUS_MPILEUP_ARGS,
    )
    await assert_md5(result.fasta, CONSENSUS_MD5, label="test.fa")
    assert result.qual.name == "test.qual.txt", result.qual.name
    assert result.mpileup is not None and result.mpileup.name == "test.mpileup"


async def run_ivar_variants(with_gff: bool, save_mpileup: bool):
    gff = await fixture("genomics/sarscov2/genome/genome.gff3") if with_gff else None
    return await ivar_variants(
        await fixture(SORTED_BAM),
        await fixture(GENOME),
        await fixture(GENOME + ".fai"),
        gff=gff,
        save_mpileup=save_mpileup,
        prefix="test",
    )


@env.task
async def test_ivar_variants_no_gff_no_mpileup() -> None:
    # upstream case: ivar/variants "test-ivar-variants-no-gff-no-mpileup"
    result = await run_ivar_variants(with_gff=False, save_mpileup=False)
    await assert_md5(result.tsv, VARIANTS_MD5, label="test.tsv")
    assert result.mpileup is None


@env.task
async def test_ivar_variants_no_gff_with_mpileup() -> None:
    # upstream case: ivar/variants "test-ivar-variants-no-gff-with-mpileup"
    result = await run_ivar_variants(with_gff=False, save_mpileup=True)
    await assert_md5(result.tsv, VARIANTS_MD5, label="test.tsv")
    assert result.mpileup is not None
    await assert_md5(result.mpileup, MPILEUP_MD5, label="test.mpileup")


@env.task
async def test_ivar_variants_with_gff_with_mpileup() -> None:
    # upstream case: ivar/variants "test-ivar-variants-with-gff-with-mpileup"
    result = await run_ivar_variants(with_gff=True, save_mpileup=True)
    await assert_md5(result.tsv, VARIANTS_GFF_MD5, label="test.tsv")
    assert result.mpileup is not None
    await assert_md5(result.mpileup, MPILEUP_MD5, label="test.mpileup")


tests = [
    test_ivar_consensus,
    test_ivar_consensus_mpileup,
    test_ivar_variants_no_gff_no_mpileup,
    test_ivar_variants_no_gff_with_mpileup,
    test_ivar_variants_with_gff_with_mpileup,
]
