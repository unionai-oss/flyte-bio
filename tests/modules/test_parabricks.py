"""Tests for flyte_bio.modules.parabricks.

The cases and expected values are upstream's parabricks tests (modules commit
9a1cc9f, Parabricks 4.8.0-1), with their test configs' ``args``. They run on
one NVIDIA GPU (a T4 on the test cluster). Upstream snapshots:

- fq2bam, fq2bammeth, rnafq2bam, applybqsr: the alignment's reads md5 (for a
  CRAM, decoded against the reference) and the index's name; fq2bammeth also
  the ``.bai``'s md5; rnafq2bam that ``Log.final.out`` reports uniquely mapped
  reads, its name, and the chimeric case's junction file name.
- haplotypecaller, deepsomatic: the output VCF's ``variantsMD5`` (md5 of the
  record lines).
- deepvariant: the first 100 records as strings (with the GVCF's name). The
  expected 100 lines are checked here by their md5 (the lines joined with
  newlines), computed from upstream's snapshot.

As upstream, the BWA index for fq2bam is built from the test FASTA by
``bwa index``, fq2bammeth's by ``bwameth.py index``, and rnafq2bam's by STAR
2.7.2a's ``genomeGenerate`` on the minigenome (with its GTF).
"""

import asyncio
import gzip
import hashlib

from flyte.io import Dir, File

from flyte_bio.modules.bwa import bwa_index
from flyte_bio.modules.bwameth import bwameth_index
from flyte_bio.modules.parabricks import (
    AlignmentResult,
    parabricks_applybqsr,
    parabricks_deepsomatic,
    parabricks_deepvariant,
    parabricks_fq2bam,
    parabricks_fq2bammeth,
    parabricks_haplotypecaller,
    parabricks_rnafq2bam,
    parabricks_star_index,
)
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_md5, env, fixture

SC2 = "genomics/sarscov2/"
HS = "genomics/homo_sapiens/"
CHR21 = HS + "genome/chr21/sequence/"
METHYLSEQ = "https://github.com/nf-core/test-datasets/raw/methylseq/"
ECOLI = "genomics/prokaryotes/escherichia_coli/"
RECALIBRATED = HS + "illumina/bam/test2.paired_end.recalibrated.sorted.bam"


async def parabricks_text(file: File) -> str:
    async with file.open("rb") as fh:
        data = bytes(await fh.read())
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode()


async def parabricks_records(file: File) -> list[str]:
    return [line for line in (await parabricks_text(file)).splitlines() if line and not line.startswith("#")]


async def parabricks_variants_md5(file: File) -> str:
    """nft-vcf's ``variantsMD5``: the md5 of a VCF's record lines, joined without newlines."""
    return hashlib.md5("".join(await parabricks_records(file)).encode()).hexdigest().lstrip("0")


async def parabricks_reads_md5(alignment: File, fasta: File | None = None) -> str:
    """nft-bam's ``getReadsMD5``: the md5 of the reads' SEQ, in file order (a CRAM decoded against ``fasta``).

    Streamed, unlike the framework's ``reads_md5``: some outputs here don't fit the test task's memory.
    """
    sam = await samtools_view(bam=alignment, fasta=fasta)
    h = hashlib.md5()
    rest = b""
    async with sam.open("rb") as fh:
        while chunk := bytes(await fh.read(1 << 20)):
            lines = (rest + chunk).split(b"\n")
            rest = lines.pop()
            for line in lines:
                h.update(line.split(b"\t")[9])
    if rest:
        h.update(rest.split(b"\t")[9])
    return h.hexdigest().lstrip("0")


async def check_reads_md5(alignment: File, expected: str, fasta: File | None = None) -> None:
    md5 = await parabricks_reads_md5(alignment, fasta)
    assert md5 == expected, f"{alignment.name} reads md5 {md5}, expected {expected}"


def check_names(result: AlignmentResult, alignment: str, index: str) -> None:
    names = (result.alignment.name, result.index.name)
    assert names == (alignment, index), f"wrote {names}"


# fq2bam


@env.task
async def test_parabricks_fq2bam_methylseq_single_end() -> None:
    # upstream case: parabricks/fq2bam "SRR389222 - fastq - se"
    fasta, reads = await asyncio.gather(
        fixture(METHYLSEQ + "reference/genome.fa"), fixture(METHYLSEQ + "testdata/SRR389222_sub1.fastq.gz")
    )
    index = await bwa_index(fasta)
    result = await parabricks_fq2bam([reads], fasta, index, single_end=True, prefix="test")
    check_names(result, "test.bam", "test.bam.bai")
    await check_reads_md5(result.alignment, "3d5b94990c7fdf90a682edb5ee0f59de")


@env.task
async def test_parabricks_fq2bam_sarscov2_paired_end() -> None:
    # upstream case: parabricks/fq2bam "sarscov2 - fastq - pe"
    fasta, r1, r2 = await asyncio.gather(
        fixture(SC2 + "genome/genome.fasta"),
        fixture(SC2 + "illumina/fastq/test_1.fastq.gz"),
        fixture(SC2 + "illumina/fastq/test_2.fastq.gz"),
    )
    index = await bwa_index(fasta)
    result = await parabricks_fq2bam([r1, r2], fasta, index, prefix="test")
    check_names(result, "test.bam", "test.bam.bai")
    await check_reads_md5(result.alignment, "2d64e4363d9f3c0e2167fce49d5087cf")


@env.task
async def test_parabricks_fq2bam_sarscov2_single_end_cram() -> None:
    # upstream case: parabricks/fq2bam "sarscov2 - fastq - se - cram"
    fasta, reads = await asyncio.gather(
        fixture(SC2 + "genome/genome.fasta"), fixture(SC2 + "illumina/fastq/test_1.fastq.gz")
    )
    index = await bwa_index(fasta)
    result = await parabricks_fq2bam([reads], fasta, index, single_end=True, output_fmt="cram", prefix="test")
    check_names(result, "test.cram", "test.cram.crai")
    await check_reads_md5(result.alignment, "30c325e1e032eb1782a280d34c0fb1c7", fasta)


# fq2bammeth


@env.task
async def test_parabricks_fq2bammeth_single_end() -> None:
    # upstream case: parabricks/fq2bammeth "SRR389222 single-end [fastq]"
    fasta, reads = await asyncio.gather(
        fixture(ECOLI + "genome/genome.fa"), fixture(ECOLI + "illumina/SRR389222_sub1.fastq.gz")
    )
    index = await bwameth_index(fasta)
    result = await parabricks_fq2bammeth([reads], fasta, index, single_end=True, prefix="test")
    check_names(result, "test.bam", "test.bam.bai")
    await check_reads_md5(result.alignment, "a7f7ca7b5eb503ab58790d64a0273ed6")
    await assert_md5(result.index, "9ec20dda17ee0441fee0598ec0e2152c", label="test.bam.bai")


# rnafq2bam

RNA_ARGS = "--num-streams-per-gpu 1 --read-files-command zcat --x3"


async def parabricks_minigenome_star_index() -> Dir:
    fasta, gtf = await asyncio.gather(fixture(HS + "genome/minigenome.fa"), fixture(HS + "genome/minigenome.gtf"))
    return await parabricks_star_index(fasta, gtf)


async def check_log_final(log: File) -> None:
    assert log.name == "test.Log.final.out", log.name
    assert "Uniquely mapped reads %" in await parabricks_text(log), "Log.final.out reports no uniquely mapped reads"


@env.task
async def test_parabricks_rnafq2bam_single_end() -> None:
    # upstream case: parabricks/rnafq2bam "homo_sapiens"
    index, fasta, reads = await asyncio.gather(
        parabricks_minigenome_star_index(),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/test_starfusion_rnaseq_1.fastq.gz"),
    )
    result = await parabricks_rnafq2bam([reads], fasta, index, single_end=True, prefix="test", args=RNA_ARGS)
    await check_log_final(result.log_final)
    assert result.bam is not None, "rna_fq2bam wrote no test.bam"
    await check_reads_md5(result.bam, "220938db5c88cdfcb9d98e10630a9e31")


@env.task
async def test_parabricks_rnafq2bam_paired_end_multiple_fastq() -> None:
    # upstream case: parabricks/rnafq2bam "homo_sapiens - paired_end - multiple_fastq"
    index, fasta, *reads = await asyncio.gather(
        parabricks_minigenome_star_index(),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/test_starfusion_rnaseq_1.fastq.gz"),
        fixture(HS + "genome/test_starfusion_rnaseq_2.fastq.gz"),
        fixture(HS + "illumina/fastq/test_rnaseq_1.fastq.gz"),
        fixture(HS + "illumina/fastq/test_rnaseq_2.fastq.gz"),
    )
    result = await parabricks_rnafq2bam(reads, fasta, index, prefix="test", args=RNA_ARGS)
    await check_log_final(result.log_final)
    assert result.bam is not None, "rna_fq2bam wrote no test.bam"
    await check_reads_md5(result.bam, "fab46f9ccf563504edc63e6a8e9cd43")


@env.task
async def test_parabricks_rnafq2bam_chimeric() -> None:
    # upstream case: parabricks/rnafq2bam "homo_sapiens_chimeric"
    index, fasta, reads = await asyncio.gather(
        parabricks_minigenome_star_index(),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/test_starfusion_rnaseq_1.fastq.gz"),
    )
    args = "--num-streams-per-gpu 1 --read-files-command zcat --out-chim-type Junctions --min-chim-segment 15 --x3"
    result = await parabricks_rnafq2bam([reads], fasta, index, single_end=True, prefix="test", args=args)
    await check_log_final(result.log_final)
    assert result.bam is not None and result.bai is not None, "rna_fq2bam wrote no test.bam/.bai"
    await check_reads_md5(result.bam, "220938db5c88cdfcb9d98e10630a9e31")
    assert result.bai.name == "test.bam.bai", result.bai.name
    assert result.junction is not None and result.junction.name == "test.Chimeric.out.junction", result.junction


# applybqsr


async def sarscov2_applybqsr_inputs() -> tuple[File, File, File, File]:
    return await asyncio.gather(
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai"),
        fixture(SC2 + "illumina/gatk/test.baserecalibrator.table"),
        fixture(SC2 + "genome/genome.fasta"),
    )


@env.task
async def test_parabricks_applybqsr_sarscov2_bam() -> None:
    # upstream case: parabricks/applybqsr "sarscov2 | paired-end | [bam]"
    bam, _, table, fasta = await sarscov2_applybqsr_inputs()
    result = await parabricks_applybqsr(bam, table, fasta, prefix="test")
    check_names(result, "test.bam", "test.bam.bai")
    await check_reads_md5(result.alignment, "894549ee3ced6b5ca2eed2563a985217")


@env.task
async def test_parabricks_applybqsr_sarscov2_cram_out() -> None:
    # upstream case: parabricks/applybqsr "sarscov2 | paired-end | [bam] | cram out"
    bam, _, table, fasta = await sarscov2_applybqsr_inputs()
    result = await parabricks_applybqsr(bam, table, fasta, output_fmt="cram", prefix="test")
    check_names(result, "test.cram", "test.cram.crai")
    await check_reads_md5(result.alignment, "894549ee3ced6b5ca2eed2563a985217", fasta)


@env.task
async def test_parabricks_applybqsr_sarscov2_intervals() -> None:
    # upstream case: parabricks/applybqsr "sarscov2 | paired-end | [bam] | intervals"
    (bam, bai, table, fasta), bed = await asyncio.gather(
        sarscov2_applybqsr_inputs(), fixture(SC2 + "genome/bed/test.bed")
    )
    result = await parabricks_applybqsr(bam, table, fasta, input_index=bai, intervals=[bed], prefix="test")
    check_names(result, "test.bam", "test.bam.bai")
    await check_reads_md5(result.alignment, "b194021b699cc5cf4b2b1f715e0b7b4c")


@env.task
async def test_parabricks_applybqsr_human_cram() -> None:
    # upstream case: parabricks/applybqsr "homo_sapiens | paired-end | [cram]"
    cram, table, fasta = await asyncio.gather(
        fixture(HS + "illumina/cram/test.paired_end.sorted.cram"),
        fixture(HS + "illumina/gatk/test.baserecalibrator.table"),
        fixture(HS + "genome/genome.fasta"),
    )
    result = await parabricks_applybqsr(cram, table, fasta, output_fmt="cram", prefix="test")
    check_names(result, "test.cram", "test.cram.crai")
    await check_reads_md5(result.alignment, "2f11e4fe3390b8ad0a1852616fd1da04", fasta)


# haplotypecaller and deepvariant, on the recalibrated chr21 BAM


async def chr21_calling_inputs(with_intervals: bool) -> tuple[File, File | None, list[File], File]:
    if not with_intervals:
        bam, fasta = await asyncio.gather(fixture(RECALIBRATED), fixture(CHR21 + "genome.fasta"))
        return bam, None, [], fasta
    bam, bai, bed, fasta = await asyncio.gather(
        fixture(RECALIBRATED),
        fixture(RECALIBRATED + ".bai"),
        fixture(CHR21 + "multi_intervals.bed"),
        fixture(CHR21 + "genome.fasta"),
    )
    return bam, bai, [bed], fasta


async def check_haplotypecaller(with_intervals: bool, args: str, name: str, md5: str) -> None:
    bam, bai, intervals, fasta = await chr21_calling_inputs(with_intervals)
    result = await parabricks_haplotypecaller(
        bam, fasta, input_index=bai, intervals=intervals, prefix="test", args=args
    )
    assert result.vcf.name == name, result.vcf.name
    variants = await parabricks_variants_md5(result.vcf)
    assert variants == md5, f"{name} variantsMD5 {variants}"


@env.task
async def test_parabricks_haplotypecaller() -> None:
    # upstream case: parabricks/haplotypecaller "human - bam"
    await check_haplotypecaller(False, "--htvc-low-memory", "test.vcf", "1cede359fcc365a936fd9665b1017cb9")


@env.task
async def test_parabricks_haplotypecaller_intervals() -> None:
    # upstream case: parabricks/haplotypecaller "human - bam - intervals"
    await check_haplotypecaller(True, "--htvc-low-memory", "test.vcf", "6e5b7197e417e4d12f6d0703b78b5caa")


@env.task
async def test_parabricks_haplotypecaller_gvcf() -> None:
    # upstream case: parabricks/haplotypecaller "human - bam - gvcf"
    await check_haplotypecaller(False, "--gvcf --htvc-low-memory", "test.g.vcf.gz", "b60554f213b9a20eb717dd010879c1b5")


@env.task
async def test_parabricks_haplotypecaller_intervals_gvcf() -> None:
    # upstream case: parabricks/haplotypecaller "human - bam - intervals - gvcf"
    await check_haplotypecaller(True, "--gvcf --htvc-low-memory", "test.g.vcf.gz", "1e36e857d1d7f5f922941075f5076b81")


async def check_deepvariant(with_intervals: bool, args: str, name: str, md5: str) -> None:
    bam, bai, intervals, fasta = await chr21_calling_inputs(with_intervals)
    result = await parabricks_deepvariant(bam, fasta, input_index=bai, intervals=intervals, prefix="test", args=args)
    assert result.vcf.name == name, result.vcf.name
    first = (await parabricks_records(result.vcf))[:100]
    actual = hashlib.md5("\n".join(first).encode()).hexdigest()
    assert actual == md5, f"{name}: first 100 records' md5 {actual}; they start {first[:3]}"


@env.task
async def test_parabricks_deepvariant() -> None:
    # upstream case: parabricks/deepvariant "human - bam"
    await check_deepvariant(False, "", "test.vcf.gz", "3a747788e5837c187ac84656b26d6bc4")


@env.task
async def test_parabricks_deepvariant_intervals() -> None:
    # upstream case: parabricks/deepvariant "human - bam - intervals"
    await check_deepvariant(True, "", "test.vcf.gz", "be48d3e7978d3a05fdb9d92175dc752a")


@env.task
async def test_parabricks_deepvariant_gvcf() -> None:
    # upstream case: parabricks/deepvariant "human - bam - gvcf"
    await check_deepvariant(False, "--gvcf", "test.g.vcf.gz", "b03fffea8d3a4fd5542c8f601fe26d25")


@env.task
async def test_parabricks_deepvariant_intervals_gvcf() -> None:
    # upstream case: parabricks/deepvariant "human - bam - intervals - gvcf"
    await check_deepvariant(True, "--gvcf", "test.g.vcf.gz", "4182ad8ed641d80028cd831f49b5ff4d")


# deepsomatic, tumour test2 against normal test


async def check_deepsomatic(prefix: str, with_intervals: bool, md5: str) -> None:
    bam = HS + "illumina/bam/"
    tumor, tumor_bai, normal, normal_bai, fasta = await asyncio.gather(
        fixture(bam + "test2.paired_end.sorted.bam"),
        fixture(bam + "test2.paired_end.sorted.bam.bai"),
        fixture(bam + "test.paired_end.sorted.bam"),
        fixture(bam + "test.paired_end.sorted.bam.bai"),
        fixture(HS + "genome/genome.fasta"),
    )
    intervals = [await fixture(HS + "genome/genome.multi_intervals.bed")] if with_intervals else []
    result = await parabricks_deepsomatic(
        tumor, normal, fasta, tumor_index=tumor_bai, normal_index=normal_bai, intervals=intervals, prefix=prefix
    )
    assert result.vcf.name == f"{prefix}.vcf.gz", result.vcf.name
    assert result.tbi is not None and result.tbi.name.endswith(".tbi"), result.tbi
    variants = await parabricks_variants_md5(result.vcf)
    assert variants == md5, f"{result.vcf.name} variantsMD5 {variants}"


@env.task
async def test_parabricks_deepsomatic() -> None:
    # upstream case: parabricks/deepsomatic "human - tumor_normal_pair"
    await check_deepsomatic("tumor_vs_normal", False, "d8cce8086db8bac0020701c15ef5311d")


@env.task
async def test_parabricks_deepsomatic_intervals() -> None:
    # upstream case: parabricks/deepsomatic "human - tumor_normal_pair - intervals"
    await check_deepsomatic("tumor_vs_normal", True, "55cfe8ffb60e1438b520ebfa16fe4486")


@env.task
async def test_parabricks_deepsomatic_gvcf() -> None:
    # upstream case: parabricks/deepsomatic "human - tumor_normal_pair - gvcf" (a prefix ending in .g, so
    # upstream's vcf and gvcf outputs are the same file)
    await check_deepsomatic("tumor_vs_normal.g", False, "d8cce8086db8bac0020701c15ef5311d")


tests = [
    test_parabricks_fq2bam_methylseq_single_end,
    test_parabricks_fq2bam_sarscov2_paired_end,
    test_parabricks_fq2bam_sarscov2_single_end_cram,
    test_parabricks_fq2bammeth_single_end,
    test_parabricks_rnafq2bam_single_end,
    test_parabricks_rnafq2bam_paired_end_multiple_fastq,
    test_parabricks_rnafq2bam_chimeric,
    test_parabricks_applybqsr_sarscov2_bam,
    test_parabricks_applybqsr_sarscov2_cram_out,
    test_parabricks_applybqsr_sarscov2_intervals,
    test_parabricks_applybqsr_human_cram,
    test_parabricks_haplotypecaller,
    test_parabricks_haplotypecaller_intervals,
    test_parabricks_haplotypecaller_gvcf,
    test_parabricks_haplotypecaller_intervals_gvcf,
    test_parabricks_deepvariant,
    test_parabricks_deepvariant_intervals,
    test_parabricks_deepvariant_gvcf,
    test_parabricks_deepvariant_intervals_gvcf,
    test_parabricks_deepsomatic,
    test_parabricks_deepsomatic_intervals,
    test_parabricks_deepsomatic_gvcf,
]
