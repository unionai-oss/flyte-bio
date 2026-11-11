"""Tests for flyte_bio.modules.longphase.

The cases and expected values are upstream's longphase/phase and
longphase/haplotag tests (modules commit 9a1cc9f, LongPhase 2.1.1), with
their test configs' ``--ont`` (phase) and ``--log --cram`` (the haplotag
"log & cram" case). The phased VCFs are compared decompressed, as nf-test
hashes them. As upstream, the haplotagged alignments are checked by their
header (whose ``@PG`` records LongPhase's command line, with upstream's 4
CPUs, the module's default) and their reads' md5.
"""

import asyncio
import hashlib

from flyte.io import File

from flyte_bio.modules.longphase import longphase_haplotag, longphase_phase
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HUMAN = "genomics/homo_sapiens/"
NANOPORE = HUMAN + "nanopore/bam/"
SNVS = HUMAN + "illumina/gvcf/test.genome.vcf.gz"
SVS = HUMAN + "illumina/vcf/NA24385_sv.vcf.gz"
READS_MD5 = "721264eb2824a3146b331f2532d10180"
INPUT_HEADER = [
    "@RG\tID:test\tSM:test",
    (
        "@PG\tID:minimap2\tPN:minimap2\tVN:2.28-r1209\tCL:minimap2 -y -x map-ont --secondary=no"
        " -R @RG\\tID:test\\tSM:test -t 30 -a genome.mmi test.bam_other.fastq.gz"
    ),
    (
        "@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 29"
        " -o test.bam_other.fastq.gz.bam --write-index"
    ),
]
HAPLOTAG_PG = "@PG\tID:longphase\tPN:longphase\tPP:samtools\tVN:2.1.1\tCL:longphase haplotag "
HAPLOTAG_CL = "--threads 4 -o test --reference genome.fasta --snp-file test.genome.vcf.gz --bam test.sorted.bam "


async def longphase_inputs(*paths: str) -> list[File]:
    return list(await asyncio.gather(*(fixture(p) for p in paths)))


async def longphase_reference() -> list[File]:
    return await longphase_inputs(HUMAN + "genome/genome.fasta", HUMAN + "genome/genome.fasta.fai")


async def longphase_header(alignment: File, fasta: File, fai: File) -> list[str]:
    header = await samtools_view(bam=alignment, fasta=fasta, fai=fai, args="-H --no-PG")
    async with header.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def longphase_reads_md5(alignment: File, fasta: File, fai: File) -> str:
    """The md5 of the reads' sequences, as nft-bam's ``getReadsMD5``."""
    sam = await samtools_view(bam=alignment, fasta=fasta, fai=fai)
    h = hashlib.md5()
    async with sam.open("rb") as fh:
        for line in bytes(await fh.read()).decode().splitlines():
            h.update(line.split("\t")[9].encode())
    return h.hexdigest().lstrip("0")


@env.task
async def test_longphase_phase() -> None:
    # upstream case: longphase/phase "[ bam, bai, snps, [], [] ], fasta, fai"
    (bam, bai, snvs), (fasta, fai) = await asyncio.gather(
        longphase_inputs(NANOPORE + "test.sorted.bam", NANOPORE + "test.sorted.bam.bai", SNVS), longphase_reference()
    )
    result = await longphase_phase([bam], [bai], snvs, fasta, fai, prefix="test", args="--ont")
    assert result.snv_vcf.path.endswith("/test.vcf.gz"), result.snv_vcf.path
    assert result.sv_vcf is None and result.mod_vcf is None, result
    await assert_gunzipped_md5(result.snv_vcf, "2abaab26f45733cb0b393c439fa9979a", label="test.vcf.gz")


@env.task
async def test_longphase_phase_svs() -> None:
    # upstream case: longphase/phase "[ bam, bai, snps, svs, [] ], fasta, fai"
    (bam, bai, snvs, svs), (fasta, fai) = await asyncio.gather(
        longphase_inputs(NANOPORE + "test.sorted.bam", NANOPORE + "test.sorted.bam.bai", SNVS, SVS),
        longphase_reference(),
    )
    result = await longphase_phase([bam], [bai], snvs, fasta, fai, prefix="test", svs=svs, args="--ont")
    assert result.sv_vcf is not None and result.mod_vcf is None, result
    await assert_gunzipped_md5(result.snv_vcf, "6714d53b06114ca758161de1072ff452", label="test.vcf.gz")
    await assert_gunzipped_md5(result.sv_vcf, "8ac6792dcc78e3ed1df18eb9b55d09ef", label="test_SV.vcf.gz")


@env.task
async def test_longphase_phase_two_bams() -> None:
    # upstream case: longphase/phase "[ bam x2, bai x2, snps, svs, [] ], fasta, fai"
    (bam, bam2, bai, bai2, snvs, svs), (fasta, fai) = await asyncio.gather(
        longphase_inputs(
            NANOPORE + "test.sorted.bam",
            NANOPORE + "test2.sorted.bam",
            NANOPORE + "test.sorted.bam.bai",
            NANOPORE + "test2.sorted.bam.bai",
            SNVS,
            SVS,
        ),
        longphase_reference(),
    )
    result = await longphase_phase([bam, bam2], [bai, bai2], snvs, fasta, fai, prefix="test", svs=svs, args="--ont")
    assert result.sv_vcf is not None, result
    await assert_gunzipped_md5(result.snv_vcf, "79701f75862c8a8b19853c41f811a291", label="test.vcf.gz")
    await assert_gunzipped_md5(result.sv_vcf, "81a8de1128e9ef6489a7cea8d00670d5", label="test_SV.vcf.gz")


@env.task
async def test_longphase_haplotag() -> None:
    # upstream case: longphase/haplotag "[ bam, bai, snps, [], [] ], fasta, fai"
    (bam, bai, snvs), (fasta, fai) = await asyncio.gather(
        longphase_inputs(NANOPORE + "test.sorted.bam", NANOPORE + "test.sorted.bam.bai", SNVS), longphase_reference()
    )
    result = await longphase_haplotag(bam, snvs, fasta, fai, prefix="test", bai=bai)
    assert result.alignment.path.endswith("/test.bam"), result.alignment.path
    assert result.log is None, result.log
    header, reads = await asyncio.gather(
        longphase_header(result.alignment, fasta, fai), longphase_reads_md5(result.alignment, fasta, fai)
    )
    expected = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr22\tLN:40001", *INPUT_HEADER, HAPLOTAG_PG + HAPLOTAG_CL]
    assert header == expected, header
    assert reads == READS_MD5, f"reads md5 {reads}"


@env.task
async def test_longphase_haplotag_log_cram() -> None:
    # upstream case: longphase/haplotag "[ bam, bai, snps, [], [] ], fasta, fai - log & cram"
    (bam, bai, snvs), (fasta, fai) = await asyncio.gather(
        longphase_inputs(NANOPORE + "test.sorted.bam", NANOPORE + "test.sorted.bam.bai", SNVS), longphase_reference()
    )
    result = await longphase_haplotag(bam, snvs, fasta, fai, prefix="test", bai=bai, args="--log --cram")
    assert result.alignment.path.endswith("/test.cram"), result.alignment.path
    assert result.log is not None, "no test.log"
    await assert_md5(result.log, "3c10627fd40507249fe963db941e99d6", label="test.log")
    header, reads = await asyncio.gather(
        longphase_header(result.alignment, fasta, fai), longphase_reads_md5(result.alignment, fasta, fai)
    )
    assert header[2:6] == [*INPUT_HEADER, HAPLOTAG_PG + "--log --cram " + HAPLOTAG_CL], header
    assert reads == READS_MD5, f"reads md5 {reads}"


@env.task
async def test_longphase_haplotag_svs() -> None:
    # upstream case: longphase/haplotag "[ bam, bai, snps, svs, [] ], fasta, fai"
    (bam, bai, snvs, svs), (fasta, fai) = await asyncio.gather(
        longphase_inputs(NANOPORE + "test.sorted.bam", NANOPORE + "test.sorted.bam.bai", SNVS, SVS),
        longphase_reference(),
    )
    result = await longphase_haplotag(bam, snvs, fasta, fai, prefix="test", bai=bai, svs=svs)
    header, reads = await asyncio.gather(
        longphase_header(result.alignment, fasta, fai), longphase_reads_md5(result.alignment, fasta, fai)
    )
    expected = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr22\tLN:40001",
        *INPUT_HEADER,
        HAPLOTAG_PG + HAPLOTAG_CL + "--sv-file NA24385_sv.vcf.gz ",
    ]
    assert header == expected, header
    assert reads == READS_MD5, f"reads md5 {reads}"


tests = [
    test_longphase_phase,
    test_longphase_phase_svs,
    test_longphase_phase_two_bams,
    test_longphase_haplotag,
    test_longphase_haplotag_log_cram,
    test_longphase_haplotag_svs,
]
