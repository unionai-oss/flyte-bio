"""Tests for flyte_bio.modules.trgt.

The cases and expected values are upstream's trgt/genotype and trgt/plot
tests (modules commit 9a1cc9f, TRGT 5.1.0). The reference is the gunzipped
chr22 genome with its samtools ``.fai``, the reads the PacBio test BAM, and
the repeats upstream's one-line BED (no trailing newline). Genotyping, with
karyotype ``XX`` and without, is checked as upstream snapshots it: the
spanning reads' sequences (nft-bam's ``getReadsMD5``) and the VCF's records
(nft-vcf's ``variantsMD5``).

The plot case runs upstream's chain: genotype → ``samtools sort`` and
``samtools index`` of the spanning reads, ``bcftools sort`` (with
``--write-index=tbi --output-type z`` and prefix ``sorted_``, its config)
of the genotypes → plot of repeat ``TEST``. As upstream, only the plot's
name is checked (``test_TEST.png``).
"""

import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bcftools import bcftools_sort
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.samtools import samtools_faidx, samtools_index, samtools_sort
from flyte_bio.modules.trgt import TrgtGenotypeResult, trgt_genotype, trgt_plot
from tests.framework import assert_reads_md5, env, fixture

BAM = "genomics/homo_sapiens/pacbio/bam/test.sorted.bam"
GENOME = "genomics/homo_sapiens/genome/chr22_chr22_KI270734v1_random/sequence/genome.fa.gz"
REPEATS_BED = "chr22\t18890357\t18890451\tID=TEST;MOTIFS=AT;STRUC=(AT)n"


async def trgt_reference() -> tuple[File, File]:
    fasta = await gunzip(await fixture(GENOME))
    fai, _ = await samtools_faidx(fasta=fasta)
    return fasta, fai


async def trgt_repeats() -> File:
    path = Path(tempfile.mkdtemp()) / "repeats.bed"
    path.write_text(REPEATS_BED)
    return await File.from_local(str(path))


async def trgt_variants_md5(vcf: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of the record lines joined without newlines, in hex without leading zeros."""
    async with vcf.open("rb") as fh:
        text = gzip.decompress(bytes(await fh.read())).decode()
    body = "".join(line for line in text.splitlines() if not line.startswith("#"))
    return format(int(hashlib.md5(body.encode()).hexdigest(), 16), "x")


async def genotype(karyotype: str) -> tuple[TrgtGenotypeResult, File, File, File]:
    fasta, fai = await trgt_reference()
    repeats = await trgt_repeats()
    result = await trgt_genotype(
        await fixture(BAM), await fixture(BAM + ".bai"), fasta, fai, repeats, karyotype=karyotype, prefix="test"
    )
    return result, fasta, fai, repeats


async def check_genotype(result: TrgtGenotypeResult) -> None:
    assert result.vcf.name == "test.vcf.gz", result.vcf.name
    assert result.spanning_bam is not None and result.spanning_bam.name == "test.spanning.bam", result.spanning_bam
    await assert_reads_md5(result.spanning_bam, "a155caa3310c034ded2fad32a6be6163", label="test.spanning.bam")
    variants = await trgt_variants_md5(result.vcf)
    assert variants == "57131f453c695cbea40812f6c246d0f", f"test.vcf.gz variantsMD5 {variants}"


@env.task
async def test_trgt_genotype_karyotype_xx() -> None:
    result, *_ = await genotype("XX")
    await check_genotype(result)


@env.task
async def test_trgt_genotype_no_karyotype() -> None:
    result, *_ = await genotype("")
    await check_genotype(result)


@env.task
async def test_trgt_plot() -> None:
    result, fasta, fai, repeats = await genotype("XX")
    assert result.spanning_bam is not None
    sorted_bam = await samtools_sort(bam=result.spanning_bam)
    bai = await samtools_index(bam=sorted_bam)
    genotypes = await bcftools_sort(result.vcf, "sorted_", args="--write-index=tbi --output-type z")
    assert genotypes.index is not None, "bcftools sort wrote no index"
    plot = await trgt_plot(
        sorted_bam, bai, genotypes.vcf, genotypes.index, "TEST", fasta, fai, repeats, prefix="test_TEST"
    )
    assert plot.name == "test_TEST.png", plot.name


tests = [
    test_trgt_genotype_karyotype_xx,
    test_trgt_genotype_no_karyotype,
    test_trgt_plot,
]
