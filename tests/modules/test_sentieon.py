"""Tests for flyte_bio.modules.sentieon — placeholders: none of them run.

Sentieon needs a licence we don't have, so every test raises
:class:`NotTested` first. The rest of each test is upstream's case (modules
commit 9a1cc9f, Sentieon 202503.03): its inputs, the call and the values
upstream's snapshot records, so deleting the ``raise`` enables it once the
tasks can reach a licence (see the module docstring for the secrets).

How upstream's snapshot values are checked once enabled:

- ``getReadsMD5`` of a BAM: :func:`tests.framework.assert_reads_md5`; of a
  CRAM: :func:`tests.modules.test_samtools.cram_reads_md5`;
- nft-vcf's ``variantsMD5``: :func:`tests.modules.test_bcftools.variants_md5`;
- nft-vcf's ``summary``: :func:`tests.modules.test_clair3.vcf_summary` (its
  chromosomes, sample count, record count and phasing);
- ``readLines()[a..b]`` slices: the same lines of the file, md5'd joined with
  newlines where upstream md5s them;
- files snapshotted by name: their names.

The Haplotyper recalibration case builds its table with
:func:`flyte_bio.modules.sentieon.sentieon_qualcal`, as upstream's setup does.
"""

import asyncio
import hashlib
import shutil
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.picard import picard_addorreplacereadgroups
from flyte_bio.modules.sentieon import (
    GermlineCallsResult,
    StarAlignResult,
    sentieon_applyvarcal,
    sentieon_bwaindex,
    sentieon_bwamem,
    sentieon_datametrics,
    sentieon_dedup,
    sentieon_dnamodelapply,
    sentieon_dnascope,
    sentieon_gvcftyper,
    sentieon_haplotyper,
    sentieon_qualcal,
    sentieon_readwriter,
    sentieon_rsemcalculateexpression,
    sentieon_rsempreparereference,
    sentieon_staralign,
    sentieon_tnscope,
    sentieon_varcal,
    sentieon_wgsmetrics,
)
from flyte_bio.modules.star import star_genome_generate
from tests.framework import NotTested, assert_md5, assert_reads_md5, env, fixture
from tests.modules.test_bcftools import variants_md5
from tests.modules.test_clair3 import vcf_summary
from tests.modules.test_samtools import cram_reads_md5

NO_LICENCE = "needs a Sentieon licence; not run"

SARS_FASTA = "genomics/sarscov2/genome/genome.fasta"
SARS_FASTQ = "genomics/sarscov2/illumina/fastq/"
SARS_BAM = "genomics/sarscov2/illumina/bam/"
HUMAN_FASTA = "genomics/homo_sapiens/genome/genome.fasta"
HUMAN_GTF = "genomics/homo_sapiens/genome/genome.gtf"
HUMAN_BED = "genomics/homo_sapiens/genome/genome.bed"
HUMAN_BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"
HUMAN_CRAM = "genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram"
HUMAN_DBSNP = "genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz"
RNASEQ_FASTQ = "genomics/homo_sapiens/illumina/fastq/test_rnaseq_"
CHR21 = "genomics/homo_sapiens/genome/chr21/"
CHR21_FASTA = CHR21 + "sequence/genome.fasta"
HAPLOTC_VCF = "genomics/homo_sapiens/illumina/gatk/haplotypecaller_calls/test2_haplotc.ann.vcf.gz"
DNASCOPE_MODEL = "https://s3.amazonaws.com/sentieon-release/other/SentieonDNAscopeModel{}.model"

# bwamem's read group, as upstream's test sets it (without its shell quotes).
READ_GROUP = r"-R @RG\tID:sample_lane\tPU:lane\tSM:patient_sample\tLB:sample\tDS:fasta\tPL:seqplatform"
# Upstream's Haplotyper test config sets these for both passes.
HAPLOTYPER_MODEL = "--genotype_model multinomial"
VARCAL_ARGS = "--var_type SNP --annotation QD --annotation MQ --annotation FS --annotation SOR"
STAR_ARGS = (
    "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outWigType bedGraph --outWigStrand Unstranded"
)
STAR_ARRIBA_ARGS = (
    "--readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0"
    " --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5"
    " --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip"
    " --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1"
    " --chimSegmentReadGapMax 3 --chimMultimapNmax 50"
)
STAR_FUSION_ARGS = (
    "--readFilesCommand zcat --outSAMtype BAM Unsorted --outReadsUnmapped None --twopassMode Basic"
    " --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8"
    " --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000"
    " --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4"
    " --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1"
    " --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30"
)


async def lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def lines_md5(file: File, first: int, last: int) -> str:
    """md5 of lines ``first..last`` (inclusive) joined with newlines, as upstream's ``readLines()[a..b].join``."""
    return hashlib.md5("\n".join((await lines(file))[first : last + 1]).encode()).hexdigest()


def assert_name(file: File | None, name: str) -> None:
    assert file is not None and file.name == name, f"expected {name}, got {file.name if file else None}"


async def check_variants(vcf: File | None, md5: str, label: str) -> None:
    assert vcf is not None, f"{label}: no VCF"
    actual = await variants_md5(vcf)
    assert actual == md5, f"{label}: variantsMD5 {actual}, expected {md5}"


async def check_summary(vcf: File | None, chromosomes: list[str], samples: int, records: int, phased: bool) -> None:
    assert vcf is not None, "no VCF"
    expected = {"chromosomes": chromosomes, "samples": samples, "records": records, "phased": phased}
    summary = await vcf_summary(vcf)
    assert summary == expected, f"vcf summary {summary} != {expected}"


async def sars_reference() -> tuple[File, File]:
    return await fixture(SARS_FASTA), await fixture(SARS_FASTA + ".fai")


async def human_reference() -> tuple[File, File]:
    return await fixture(HUMAN_FASTA), await fixture(HUMAN_FASTA + ".fai")


async def chr21_reference() -> tuple[File, File]:
    return await fixture(CHR21_FASTA), await fixture(CHR21_FASTA + ".fai")


# --- bwaindex ---


@env.task
async def test_sentieon_bwaindex_fasta() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/bwaindex "build index - fasta"
    index = await sentieon_bwaindex(await fixture(SARS_FASTA))
    expected = {
        "genome.amb": "3a68b8b2287e07dd3f5f95f4344ba76e",
        "genome.ann": "c32e11f6c859f166c7525a9c1d583567",
        "genome.bwt": "0469c30a1e239dd08f68afe66fde99da",
        "genome.pac": "983e3d2cd6f36e2546e6d25a0da78d66",
        "genome.sa": "ab3952cabf026b48cd3eb5bccbb636d1",
    }
    for name, md5 in expected.items():
        found = await index.get_file(name)
        assert found is not None, f"no {name}"
        await assert_md5(found, md5, label=name)


# --- bwamem ---


async def bwamem_case(paired: bool, prefix: str, reads_md5: str, index_name: str) -> None:
    fasta, fai = await sars_reference()
    index = await sentieon_bwaindex(fasta)
    reads_2 = await fixture(SARS_FASTQ + "test_2.fastq.gz") if paired else None
    out = await sentieon_bwamem(
        await fixture(SARS_FASTQ + "test_1.fastq.gz"), index, fasta, prefix, reads_2=reads_2, fai=fai, args=READ_GROUP
    )
    if prefix.endswith(".cram"):
        actual = await cram_reads_md5(out.alignment, fasta, fai)
        assert actual == reads_md5, f"reads md5 {actual}, expected {reads_md5}"
    else:
        await assert_reads_md5(out.alignment, reads_md5, label="bwa mem")
    assert_name(out.index, index_name)


@env.task
async def test_sentieon_bwamem_single_end() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/bwamem "Single-End"
    await bwamem_case(False, "test.bam", "30c325e1e032eb1782a280d34c0fb1c7", "test.bam.bai")


@env.task
async def test_sentieon_bwamem_single_end_cram() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/bwamem "Single-End Output CRAM"
    await bwamem_case(False, "test.cram", "30c325e1e032eb1782a280d34c0fb1c7", "test.cram.crai")


@env.task
async def test_sentieon_bwamem_paired_end() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/bwamem "Paired-End"
    await bwamem_case(True, "test.bam", "2d64e4363d9f3c0e2167fce49d5087cf", "test.bam.bai")


# --- readwriter ---

SARS_MERGE_BAMS = [
    SARS_BAM + "test.paired_end.methylated.sorted.bam",
    SARS_BAM + "test.paired_end.sorted.bam",
    SARS_BAM + "test.single_end.sorted.bam",
]


async def readwriter_case(
    paths: list[str], index_suffix: str, prefix: str, reference: tuple[File, File] | None, reads_md5: str
) -> None:
    alignments = await asyncio.gather(*(fixture(p) for p in paths))
    indexes = await asyncio.gather(*(fixture(p + index_suffix) for p in paths))
    fasta, fai = reference or (None, None)
    out = await sentieon_readwriter(list(alignments), list(indexes), prefix, fasta=fasta, fai=fai)
    if prefix.endswith(".cram"):
        assert fasta is not None and fai is not None
        actual = await cram_reads_md5(out.alignment, fasta, fai)
        assert actual == reads_md5, f"reads md5 {actual}, expected {reads_md5}"
    else:
        await assert_reads_md5(out.alignment, reads_md5, label="ReadWriter")
    assert_name(out.index, f"{prefix}.bai" if prefix.endswith(".bam") else f"{prefix}.crai")


@env.task
async def test_sentieon_readwriter_merge_bam() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/readwriter "merge bam" (one BAM, no reference)
    await readwriter_case([HUMAN_BAM], ".bai", "test.bam", None, "2f11e4fe3390b8ad0a1852616fd1da04")


@env.task
async def test_sentieon_readwriter_cram() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/readwriter "cram"
    await readwriter_case(
        [HUMAN_CRAM], ".crai", "test.cram", await human_reference(), "2f11e4fe3390b8ad0a1852616fd1da04"
    )


@env.task
async def test_sentieon_readwriter_merging_bam() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/readwriter "merging - output bam"
    await readwriter_case(SARS_MERGE_BAMS, ".bai", "test.bam", None, "ed89d81d2d6adafdbb819a1a4ec86493")


@env.task
async def test_sentieon_readwriter_merging_bam_to_cram() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/readwriter "merging - output bam to cram"
    await readwriter_case(
        SARS_MERGE_BAMS, ".bai", "test.cram", await sars_reference(), "ed89d81d2d6adafdbb819a1a4ec86493"
    )


# --- dedup ---


async def dedup_case(prefix: str, args4: str) -> None:
    fasta, fai = await sars_reference()
    bam = await fixture(SARS_BAM + "test.paired_end.sorted.bam")
    bai = await fixture(SARS_BAM + "test.paired_end.sorted.bam.bai")
    out = await sentieon_dedup([bam], [bai], fasta, fai, prefix, args4=args4)
    if prefix.endswith(".cram"):
        actual = await cram_reads_md5(out.alignment, fasta, fai)
        assert actual == "894549ee3ced6b5ca2eed2563a985217", f"reads md5 {actual}"
        assert_name(out.index, f"{prefix}.crai")
    else:
        await assert_reads_md5(out.alignment, "894549ee3ced6b5ca2eed2563a985217", label="Dedup")
        assert_name(out.index, f"{prefix}.bai")
    # Upstream snapshots that none of these mention SentieonCommandLine.
    for file in (out.score, out.metrics, out.metrics_multiqc_tsv):
        assert not any("SentieonCommandLine" in line for line in await lines(file)), f"{file.name} names the command"


@env.task
async def test_sentieon_dedup_mark_bam() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/dedup "Test marking duplicates, bam"
    await dedup_case("test.dedup.bam", "")


@env.task
async def test_sentieon_dedup_remove_cram() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/dedup "Test removing duplicates, cram"
    await dedup_case("test.dedup.cram", "--rmdup")


# --- datametrics ---


async def datametrics_case(plot_results: bool) -> None:
    fasta, fai = await sars_reference()
    bam = await fixture(SARS_BAM + "test.paired_end.sorted.bam")
    bai = await fixture(SARS_BAM + "test.paired_end.sorted.bam.bai")
    out = await sentieon_datametrics([bam], [bai], fasta, fai, "test", plot_results=plot_results)
    expected = [
        (out.mq_metrics, "8fdd2f19b222aa3a8f1d4f82d3792809"),
        (out.qd_metrics, "a5b70394b8ab8ab26fad4a0b819df96c"),
        (out.gc_summary, "1e40f4ebb6098e8f7e07f51d873cc26b"),
        (out.gc_metrics, "547b7f29acc4efc43462d4df3f34b774"),
        (out.aln_metrics, "8d28671f4ef0d56ed1704cf236bf0952"),
        (out.is_metrics, "186761dcce9ac40f25ae4bdaa0153544"),
    ]
    for file, md5 in expected:
        actual = await lines_md5(file, 1, 3)
        assert actual == md5, f"{file.name} lines 2-4: md5 {actual}, expected {md5}"
    if plot_results:
        assert_name(out.mq_plot, "test_mq_metrics.pdf")
        assert_name(out.qd_plot, "test_qd_metrics.pdf")
        assert_name(out.is_plot, "test_is_metrics.pdf")
        assert_name(out.gc_plot, "test_gc_metrics.pdf")
    else:
        assert out.mq_plot is None and out.qd_plot is None and out.is_plot is None and out.gc_plot is None


@env.task
async def test_sentieon_datametrics_no_plots() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/datametrics "metrics - no plots"
    await datametrics_case(False)


@env.task
async def test_sentieon_datametrics_with_plots() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/datametrics "metrics - with plots"
    await datametrics_case(True)


# --- wgsmetrics ---


@env.task
async def test_sentieon_wgsmetrics_bam() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/wgsmetrics "sarscov2 - bam"
    fasta, fai = await sars_reference()
    bam = await fixture(SARS_BAM + "test.paired_end.sorted.bam")
    bai = await fixture(SARS_BAM + "test.paired_end.sorted.bam.bai")
    out = await sentieon_wgsmetrics([bam], [bai], fasta, fai, "test")
    expected = [
        (
            "GENOME_TERRITORY\tMEAN_COVERAGE\tSD_COVERAGE\tMEDIAN_COVERAGE\tMAD_COVERAGE\tPCT_EXC_MAPQ\tPCT_EXC_DUPE\t"
            "PCT_EXC_UNPAIRED\tPCT_EXC_BASEQ\tPCT_EXC_OVERLAP\tPCT_EXC_CAPPED\tPCT_EXC_TOTAL\tPCT_1X\tPCT_5X\tPCT_10X\t"
            "PCT_15X\tPCT_20X\tPCT_25X\tPCT_30X\tPCT_40X\tPCT_50X\tPCT_60X\tPCT_70X\tPCT_80X\tPCT_90X\tPCT_100X\t"
            "HET_SNP_SENSITIVITY\tHET_SNP_Q\t"
        ),
        (
            "29829\t0.63086\t0.94268\t0\t0\t0.0054378\t0\t0.013649\t0.063246\t0.2309\t0\t0.31324\t0.41151\t0.0051963\t"
            "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0.17212\t1\t"
        ),
        "",
        "coverage\thigh_quality_coverage_count\t",
        "0\t17554",
    ]
    actual = (await lines(out))[1:6]
    assert actual == expected, f"wgs metrics lines 2-6: {actual}"


# --- haplotyper ---


async def haplotyper_case(
    emit_vcf: str,
    emit_gvcf: bool,
    vcf_md5: str | None,
    gvcf_md5: str | None,
    intervals: bool = False,
    dbsnp: bool = False,
    recal_table: File | None = None,
    alignments: tuple[list[str], list[str]] | None = None,
) -> None:
    fasta, fai = await human_reference()
    paths, index_paths = alignments or ([HUMAN_CRAM], [HUMAN_CRAM + ".crai"])
    out: GermlineCallsResult = await sentieon_haplotyper(
        list(await asyncio.gather(*(fixture(p) for p in paths))),
        list(await asyncio.gather(*(fixture(p) for p in index_paths))),
        fasta,
        fai,
        "test",
        emit_vcf=emit_vcf,
        emit_gvcf=emit_gvcf,
        intervals=await fixture(HUMAN_BED) if intervals else None,
        recal_table=recal_table,
        dbsnp=await fixture(HUMAN_DBSNP) if dbsnp else None,
        dbsnp_tbi=await fixture(HUMAN_DBSNP + ".tbi") if dbsnp else None,
        args2=HAPLOTYPER_MODEL,
        args3=HAPLOTYPER_MODEL,
    )
    if gvcf_md5 is not None:
        assert_name(out.gvcf_tbi, "test.g.vcf.gz.tbi")
        await check_variants(out.gvcf, gvcf_md5, "gvcf")
    else:
        assert out.gvcf is None
    if vcf_md5 is not None:
        assert_name(out.vcf_tbi, "test.unfiltered.vcf.gz.tbi")
        await check_variants(out.vcf, vcf_md5, "vcf")
    else:
        assert out.vcf is None


@env.task
async def test_sentieon_haplotyper_vcf() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/haplotyper "Sentieon Haplotyper VCF"
    await haplotyper_case("variant", False, "cea0045051da7877b38a1e25df812a91", None)


@env.task
async def test_sentieon_haplotyper_gvcf() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/haplotyper "Sentieon Haplotyper GVCF"
    await haplotyper_case("", True, None, "338fc3c37b208d6595948576833eb665")


@env.task
async def test_sentieon_haplotyper_both() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/haplotyper "Sentieon Haplotyper BOTH"
    await haplotyper_case("variant", True, "cea0045051da7877b38a1e25df812a91", "338fc3c37b208d6595948576833eb665")


@env.task
async def test_sentieon_haplotyper_intervals_both() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/haplotyper "Sentieon Haplotyper Intervals BOTH"
    await haplotyper_case(
        "variant", True, "cea0045051da7877b38a1e25df812a91", "338fc3c37b208d6595948576833eb665", intervals=True
    )


@env.task
async def test_sentieon_haplotyper_dbsnp_both() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/haplotyper "Sentieon Haplotyper DBSNP BOTH"
    await haplotyper_case(
        "variant",
        True,
        "cc1f3d4bd615f3640e7fd103cc39d2f8",
        "228556b7921205f023fec51098feeb97",
        intervals=True,
        dbsnp=True,
    )


@env.task
async def test_sentieon_haplotyper_recalibration() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/haplotyper "Sentieon Haplotyper Recalibration"; its setup's QualCal table
    fasta, fai = await human_reference()
    qualcal = await sentieon_qualcal(
        [await fixture(HUMAN_CRAM)],
        [await fixture(HUMAN_CRAM + ".crai")],
        fasta,
        fai,
        "test",
        known_sites=[await fixture(HUMAN_DBSNP)],
        known_sites_tbi=[await fixture(HUMAN_DBSNP + ".tbi")],
    )
    assert qualcal.table is not None, "QualCal wrote no table"
    await haplotyper_case("variant", False, "10faa3b669c49826098e09784d8a4716", None, recal_table=qualcal.table)


@env.task
async def test_sentieon_haplotyper_multiple_crams() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/haplotyper "Sentieon Haplotyper multiple CRAMs"
    bam2 = "genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam"
    await haplotyper_case(
        "variant",
        False,
        "b5d6e09e336438e38f7bf5531799e3a",
        None,
        alignments=([HUMAN_CRAM, bam2], [HUMAN_CRAM + ".crai", bam2 + ".bai"]),
    )


# --- gvcftyper ---

GVCF = "genomics/homo_sapiens/illumina/gvcf/test.genome.vcf"


async def gvcftyper_case(gvcf: str, index: str, md5: str, intervals: bool = False, dbsnp: bool = False) -> None:
    fasta, fai = await human_reference()
    out = await sentieon_gvcftyper(
        [await fixture(gvcf)],
        [await fixture(index)],
        fasta,
        fai,
        "test_genotyped",
        intervals=await fixture(HUMAN_BED) if intervals else None,
        dbsnp=await fixture(HUMAN_DBSNP) if dbsnp else None,
        dbsnp_tbi=await fixture(HUMAN_DBSNP + ".tbi") if dbsnp else None,
    )
    assert_name(out.tbi, "test_genotyped.vcf.gz.tbi")
    await check_variants(out.vcf, md5, "GVCFtyper")


@env.task
async def test_sentieon_gvcftyper_vcf() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/gvcftyper "sentieon gvcftyper vcf"
    await gvcftyper_case(GVCF, GVCF + ".idx", "d13216836f1452e200b215b796606671")


@env.task
async def test_sentieon_gvcftyper_vcf_gz() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/gvcftyper "sentieon gvcftyper vcf.gz"
    await gvcftyper_case(GVCF + ".gz", GVCF + ".gz.tbi", "d13216836f1452e200b215b796606671")


@env.task
async def test_sentieon_gvcftyper_dbsnp() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/gvcftyper "sentieon gvcftyper dbsnp"
    await gvcftyper_case(GVCF + ".gz", GVCF + ".gz.tbi", "21606383c760bf676d4c1f747b97d118", dbsnp=True)


@env.task
async def test_sentieon_gvcftyper_intervals() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/gvcftyper "sentieon gvcftyper intervals"
    await gvcftyper_case(GVCF + ".gz", GVCF + ".gz.tbi", "d13216836f1452e200b215b796606671", intervals=True)


@env.task
async def test_sentieon_gvcftyper_dbsnp_intervals() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/gvcftyper "sentieon gvcftyper dbsnp intervals"
    await gvcftyper_case(GVCF + ".gz", GVCF + ".gz.tbi", "21606383c760bf676d4c1f747b97d118", intervals=True, dbsnp=True)


# --- dnascope ---


async def dnascope_calls(model_version: str | None, pcr_indel_model: str, emit_vcf: str, emit_gvcf: bool):
    fasta, fai = await human_reference()
    return await sentieon_dnascope(
        await fixture(HUMAN_BAM),
        await fixture(HUMAN_BAM + ".bai"),
        fasta,
        fai,
        "test",
        emit_vcf=emit_vcf,
        emit_gvcf=emit_gvcf,
        pcr_indel_model=pcr_indel_model,
        ml_model=await fixture(DNASCOPE_MODEL.format(model_version)) if model_version else None,
    )


@env.task
async def test_sentieon_dnascope_conservative_variant() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/dnascope "conservative indel model - variant"
    out = await dnascope_calls("1.0", "CONSERVATIVE", "variant", False)
    await check_summary(out.vcf, ["chr22"], 1, 5, False)
    assert_name(out.vcf_tbi, "test.unfiltered.vcf.gz.tbi")


@env.task
async def test_sentieon_dnascope_conservative_all() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/dnascope "conservative indel model - all" (upstream passes AGGRESSIVE)
    out = await dnascope_calls("1.0", "AGGRESSIVE", "all", False)
    await check_summary(out.vcf, ["chr22"], 1, 8, False)
    assert_name(out.vcf_tbi, "test.unfiltered.vcf.gz.tbi")


@env.task
async def test_sentieon_dnascope_hostile_confident() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/dnascope "hostile indel model - confident" (no model; a GVCF too)
    out = await dnascope_calls(None, "HOSTILE", "confident", True)
    await check_summary(out.vcf, ["chr22"], 1, 5, False)
    assert_name(out.vcf_tbi, "test.unfiltered.vcf.gz.tbi")


# --- dnamodelapply ---


@env.task
async def test_sentieon_dnamodelapply_apply_model() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/dnamodelapply "apply model", on its setup's DNAscope calls
    calls = await dnascope_calls("1.1", "CONSERVATIVE", "variant", False)
    assert calls.vcf is not None and calls.vcf_tbi is not None
    fasta, fai = await human_reference()
    model = await fixture(DNASCOPE_MODEL.format("1.1"))
    out = await sentieon_dnamodelapply(calls.vcf, calls.vcf_tbi, fasta, fai, model, "test_applied")
    await check_summary(out.vcf, ["chr22"], 1, 5, False)
    assert_name(out.tbi, "test_applied.vcf.gz.tbi")


# --- tnscope ---


async def tnscope_case(dbsnp: bool) -> None:
    recalibrated = "genomics/homo_sapiens/illumina/bam/{}.paired_end.recalibrated.sorted.bam"
    normal = await picard_addorreplacereadgroups(
        await fixture(recalibrated.format("test")),
        "normal",
        args="--CREATE_INDEX true -ID normal -LB normal -PL ILLUMINA -PU bc1 -SM normal",
    )
    tumour = await picard_addorreplacereadgroups(
        await fixture(recalibrated.format("test2")),
        "tumour",
        args="--CREATE_INDEX true -ID tumour -LB tumour -PL ILLUMINA -PU bc1 -SM tumour",
    )
    assert normal.bai is not None and tumour.bai is not None
    fasta, fai = await chr21_reference()
    germline = CHR21 + "germlineresources/dbsnp_138.hg38.vcf.gz"
    out = await sentieon_tnscope(
        [normal.reads, tumour.reads],
        [normal.bai, tumour.bai],
        fasta,
        fai,
        "tumour",
        dbsnp=await fixture(germline) if dbsnp else None,
        dbsnp_tbi=await fixture(germline + ".tbi") if dbsnp else None,
        args2="--tumor_sample tumour --normal_sample normal",
    )
    # Upstream's summary: no chromosomes, 2 samples, 0 variants. nft-vcf reports phased=true for a VCF
    # without records; a record-less VCF has no phased genotype to find, hence False here.
    await check_summary(out.vcf, [], 2, 0, False)
    assert_name(out.tbi, "tumour.vcf.gz.tbi")


@env.task
async def test_sentieon_tnscope_tumour_normal() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/tnscope "tumour normal - no dbsnp"
    await tnscope_case(False)


@env.task
async def test_sentieon_tnscope_tumour_normal_dbsnp() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/tnscope "tumour normal - dbsnp"
    await tnscope_case(True)


# --- varcal ---

GERMLINE = CHR21 + "germlineresources/"
RECAL_HEADER = [
    "##fileformat=VCFv4.2",
    '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">',
    (
        '##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the negative'
        ' training set of bad variants">'
    ),
    (
        '##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the positive'
        ' training set of good variants">'
    ),
    (
        '##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false'
        ' under the trained gaussian mixture model">'
    ),
    (
        '##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the'
        ' Gaussian mixture model, likely the reason why the variant was filtered out">'
    ),
]
HAPMAP_LABEL = "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz"
OMNI_LABEL = "--resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz"
PHASE1_LABEL = "--resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.hg38.vcf.gz"
DBSNP_LABEL = "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf.gz"


async def varcal_run(resources: list[str], labels: list[str]):
    fasta, fai = await chr21_reference()
    return await sentieon_varcal(
        await fixture(HAPLOTC_VCF),
        await fixture(HAPLOTC_VCF + ".tbi"),
        list(await asyncio.gather(*(fixture(GERMLINE + r) for r in resources))),
        list(await asyncio.gather(*(fixture(GERMLINE + r + ".tbi") for r in resources))),
        labels,
        fasta,
        fai,
        "test",
        args=VARCAL_ARGS,
    )


async def varcal_case(resources: list[str], labels: list[str], tranches_md5: str) -> None:
    out = await varcal_run(resources, labels)
    header = (await lines(out.recal))[0:6]
    assert header == RECAL_HEADER, f"recal header {header}"
    assert_name(out.idx, "test.recal.idx")
    await assert_md5(out.tranches, tranches_md5, label="tranches")


@env.task
async def test_sentieon_varcal_one_resource() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/varcal "homo sapiens - one resource"
    await varcal_case(["hapmap_3.3.hg38.vcf.gz"], [HAPMAP_LABEL], "21e65cd1b310b4f51fa667856a59a17a")


@env.task
async def test_sentieon_varcal_multiple_resources() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/varcal "homo sapiens - multiple resources"
    await varcal_case(
        [
            "hapmap_3.3.hg38.vcf.gz",
            "1000G_omni2.5.hg38.vcf.gz",
            "1000G_phase1.snps.hg38.vcf.gz",
            "dbsnp_138.hg38.vcf.gz",
        ],
        [HAPMAP_LABEL, OMNI_LABEL, PHASE1_LABEL, DBSNP_LABEL],
        "575517ff6abbc8d7bf7aee4eef449d29",
    )


@env.task
async def test_sentieon_varcal_packed_resources() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/varcal "homo sapiens - list element packs multiple --resource directives";
    # upstream checks only that it succeeds
    await varcal_run(
        ["hapmap_3.3.hg38.vcf.gz", "1000G_omni2.5.hg38.vcf.gz", "dbsnp_138.hg38.vcf.gz"],
        [f"{HAPMAP_LABEL} {OMNI_LABEL}", DBSNP_LABEL],
    )


# --- applyvarcal ---


@env.task
async def test_sentieon_applyvarcal_homo_sapiens() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/applyvarcal "homo sapiens"
    recal = "genomics/homo_sapiens/illumina/gatk/variantrecalibrator/test2"
    fasta, fai = await chr21_reference()
    out = await sentieon_applyvarcal(
        await fixture(HAPLOTC_VCF),
        await fixture(HAPLOTC_VCF + ".tbi"),
        await fixture(recal + ".recal"),
        await fixture(recal + ".recal.idx"),
        await fixture(recal + ".tranches"),
        fasta,
        fai,
        "test_applyvarcal",
        args2="--var_type SNP --sensitivity 99.0",
    )
    await check_summary(out.vcf, ["chr21"], 1, 733, False)
    assert_name(out.tbi, "test_applyvarcal.vcf.gz.tbi")


# --- staralign ---


async def staralign_run(reads_1: list[File], reads_2: list[File], args: str) -> StarAlignResult:
    index = await star_genome_generate(
        fasta=await fixture(HUMAN_FASTA), gtf=await fixture(HUMAN_GTF), args="--genomeSAindexNbases 9"
    )
    return await sentieon_staralign(reads_1, index, "test", reads_2=reads_2, gtf=await fixture(HUMAN_GTF), args=args)


async def check_star_common(out: StarAlignResult, sj_md5: str) -> None:
    assert_name(out.log_final, "test.Log.final.out")
    assert_name(out.log_out, "test.Log.out")
    assert_name(out.log_progress, "test.Log.progress.out")
    assert out.fastq == [] and out.read_per_gene_tab is None and out.sam is None and out.wig == []
    assert out.spl_junc_tab is not None
    await assert_md5(out.spl_junc_tab, sj_md5, label="SJ.out.tab")
    assert [f.name for f in out.tab] == ["test.SJ.out.tab"], [f.name for f in out.tab]


async def check_star_signal(out: StarAlignResult, reads_md5: str, unique_md5: str, multiple_md5: str) -> None:
    assert len(out.bam) == 1 and out.bam_sorted_aligned is not None, [f.name for f in out.bam]
    await assert_reads_md5(out.bam[0], reads_md5, label="bam")
    await assert_reads_md5(out.bam_sorted_aligned, reads_md5, label="sorted bam")
    names = [f.name for f in out.bedgraph]
    assert names == ["test.Signal.Unique.str1.out.bg", "test.Signal.UniqueMultiple.str1.out.bg"], names
    await assert_md5(out.bedgraph[0], unique_md5, label="Signal.Unique")
    await assert_md5(out.bedgraph[1], multiple_md5, label="Signal.UniqueMultiple")


@env.task
async def test_sentieon_staralign_single_end() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/staralign "homo_sapiens - single_end"
    out = await staralign_run([await fixture(RNASEQ_FASTQ + "1.fastq.gz")], [], STAR_ARGS)
    await check_star_signal(
        out, "9f76be49a6607613a64f760101bdddce", "c56fc1472776fb927eaf62d973da5f9a", "e93373cf6f2a2a9506e2efdb260cdd4f"
    )
    await check_star_common(out, "75a516ab950fb958f40b29996474949c")


@env.task
async def test_sentieon_staralign_paired_end() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/staralign "homo_sapiens - paired_end"
    out = await staralign_run(
        [await fixture(RNASEQ_FASTQ + "1.fastq.gz")], [await fixture(RNASEQ_FASTQ + "2.fastq.gz")], STAR_ARGS
    )
    await check_star_signal(
        out, "db9a8324b5163b025bcc0c33e848486", "d7bf8b70b436ca048a62513e1d0ece3a", "686d58493b9eb445b56ace4d67f76ef6"
    )
    await check_star_common(out, "844af19ab0fc8cd9a3f75228445aca0d")


@env.task
async def test_sentieon_staralign_paired_end_arriba() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/staralign "homo_sapiens - paired_end - arriba"
    out = await staralign_run(
        [await fixture(RNASEQ_FASTQ + "1.fastq.gz")], [await fixture(RNASEQ_FASTQ + "2.fastq.gz")], STAR_ARRIBA_ARGS
    )
    assert len(out.bam) == 1, [f.name for f in out.bam]
    await assert_reads_md5(out.bam[0], "1a3abe88fb2490589c58497d39921bcc", label="bam")
    assert out.bedgraph == []
    await check_star_common(out, "5155c9fd1f787ad6d7d80987fb06219c")


@env.task
async def test_sentieon_staralign_paired_end_starfusion() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/staralign "homo_sapiens - paired_end - starfusion"
    out = await staralign_run(
        [await fixture(RNASEQ_FASTQ + "1.fastq.gz")], [await fixture(RNASEQ_FASTQ + "2.fastq.gz")], STAR_FUSION_ARGS
    )
    assert_name(out.junction, "test.Chimeric.out.junction")
    assert len(out.bam) == 1, [f.name for f in out.bam]
    await assert_reads_md5(out.bam[0], "caee9dcda13882d4913456973c25b57a", label="bam")
    assert out.bedgraph == []
    await check_star_common(out, "19c3faa1bfa9a0cc5e4c45f17065b53a")


async def copy_as(local: str, name: str) -> File:
    """A copy of ``local`` named ``name``, so two copies of one fixture can be staged side by side."""
    target = Path(tempfile.mkdtemp()) / name
    shutil.copyfile(local, target)
    return await File.from_local(str(target))


@env.task
async def test_sentieon_staralign_paired_end_multiple() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/staralign "homo_sapiens - paired_end - multiple": each mate's fixture twice.
    # A list's files are staged under their names, so the second copies are renamed (they sort first,
    # which doesn't change the reads STAR sees).
    r1 = await fixture(RNASEQ_FASTQ + "1.fastq.gz")
    r2 = await fixture(RNASEQ_FASTQ + "2.fastq.gz")
    r1_copy = await copy_as(await r1.download(), "test_rnaseq_1.again.fastq.gz")
    r2_copy = await copy_as(await r2.download(), "test_rnaseq_2.again.fastq.gz")
    out = await staralign_run([r1, r1_copy], [r2, r2_copy], STAR_ARGS)
    await check_star_signal(
        out, "3e54e45f5dc3e9c1f2fc55bc41531a87", "d7bf8b70b436ca048a62513e1d0ece3a", "686d58493b9eb445b56ace4d67f76ef6"
    )
    await check_star_common(out, "069877e053714e23010fe4e1c003b4a2")


# --- rsempreparereference ---


@env.task
async def test_sentieon_rsempreparereference_homo_sapiens() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/rsempreparereference "homo_sapiens" (RSEM only: the test sets no --star)
    out = await sentieon_rsempreparereference(await fixture(HUMAN_FASTA), await fixture(HUMAN_GTF))
    expected = {
        "genome.chrlist": "b190587cae0531f3cf25552d8aa674db",
        "genome.fasta": "f315020d899597c1b57e5fe9f60f4c3e",
        "genome.grp": "c2848a8b6d495956c11ec53efc1de67e",
        "genome.idx.fa": "050c521a2719c2ae48267c1e65218f29",
        "genome.n2g.idx.fa": "050c521a2719c2ae48267c1e65218f29",
        "genome.seq": "94da0c6b88c33e63c9a052a11f4f57c1",
        "genome.ti": "c9e4ae8d4d13a504eec2acf1b8589a66",
        "genome.transcripts.fa": "050c521a2719c2ae48267c1e65218f29",
    }
    names = sorted([f.name async for f in out.index.walk()])
    assert names == sorted(expected), f"index holds {names}"
    for name, md5 in expected.items():
        found = await out.index.get_file(name)
        assert found is not None, f"no {name}"
        await assert_md5(found, md5, label=name)
    await assert_md5(out.transcript_fasta, "050c521a2719c2ae48267c1e65218f29", label="transcripts")


# --- rsemcalculateexpression ---


async def rsemcalculateexpression_case(reads_1: File, reads_2: File | None) -> None:
    reference = await sentieon_rsempreparereference(await fixture(HUMAN_FASTA), await fixture(HUMAN_GTF), args="--star")
    out = await sentieon_rsemcalculateexpression(
        reads_1,
        reference.index,
        "test",
        reads_2=reads_2,
        strandedness="forward",
        args="--star --star-gzipped-read-file",
    )
    await assert_md5(out.counts_gene, "c7ec226f76736ea805771e73553ae359", label="genes.results")
    await assert_md5(out.counts_transcript, "99f7f80aa505b44ca429fdebbd7dd5d8", label="isoforms.results")
    expected = {
        "test.cnt": "76249e6b2f3c104f414aae596ba2c2f4",
        "test.model": "a7a4bc1734918ef5848604e3362b83e2",
        "test.theta": "de2e4490c98cc5383a86ae8225fd0a28",
    }
    assert [f.name for f in out.stat] == sorted(expected), [f.name for f in out.stat]
    for file in out.stat:
        await assert_md5(file, expected[file.name], label=file.name)


@env.task
async def test_sentieon_rsemcalculateexpression_bam() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/rsemcalculateexpression "homo_sapiens - bam"
    await rsemcalculateexpression_case(await fixture("genomics/homo_sapiens/illumina/bam/rsem.transcript.bam"), None)


@env.task
async def test_sentieon_rsemcalculateexpression_fastq() -> None:
    raise NotTested(NO_LICENCE)
    # upstream case: sentieon/rsemcalculateexpression "homo_sapiens - fastq" (upstream also checks a log exists)
    await rsemcalculateexpression_case(
        await fixture(RNASEQ_FASTQ + "1.fastq.gz"), await fixture(RNASEQ_FASTQ + "2.fastq.gz")
    )


tests = [
    test_sentieon_bwaindex_fasta,
    test_sentieon_bwamem_single_end,
    test_sentieon_bwamem_single_end_cram,
    test_sentieon_bwamem_paired_end,
    test_sentieon_readwriter_merge_bam,
    test_sentieon_readwriter_cram,
    test_sentieon_readwriter_merging_bam,
    test_sentieon_readwriter_merging_bam_to_cram,
    test_sentieon_dedup_mark_bam,
    test_sentieon_dedup_remove_cram,
    test_sentieon_datametrics_no_plots,
    test_sentieon_datametrics_with_plots,
    test_sentieon_wgsmetrics_bam,
    test_sentieon_haplotyper_vcf,
    test_sentieon_haplotyper_gvcf,
    test_sentieon_haplotyper_both,
    test_sentieon_haplotyper_intervals_both,
    test_sentieon_haplotyper_dbsnp_both,
    test_sentieon_haplotyper_recalibration,
    test_sentieon_haplotyper_multiple_crams,
    test_sentieon_gvcftyper_vcf,
    test_sentieon_gvcftyper_vcf_gz,
    test_sentieon_gvcftyper_dbsnp,
    test_sentieon_gvcftyper_intervals,
    test_sentieon_gvcftyper_dbsnp_intervals,
    test_sentieon_dnascope_conservative_variant,
    test_sentieon_dnascope_conservative_all,
    test_sentieon_dnascope_hostile_confident,
    test_sentieon_dnamodelapply_apply_model,
    test_sentieon_tnscope_tumour_normal,
    test_sentieon_tnscope_tumour_normal_dbsnp,
    test_sentieon_varcal_one_resource,
    test_sentieon_varcal_multiple_resources,
    test_sentieon_varcal_packed_resources,
    test_sentieon_applyvarcal_homo_sapiens,
    test_sentieon_staralign_single_end,
    test_sentieon_staralign_paired_end,
    test_sentieon_staralign_paired_end_arriba,
    test_sentieon_staralign_paired_end_starfusion,
    test_sentieon_staralign_paired_end_multiple,
    test_sentieon_rsempreparereference_homo_sapiens,
    test_sentieon_rsemcalculateexpression_bam,
    test_sentieon_rsemcalculateexpression_fastq,
]
