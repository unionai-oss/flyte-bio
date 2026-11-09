"""Tests for flyte_bio.modules.strelka.

The cases are upstream's strelka/somatic tests (modules commit 9a1cc9f,
Strelka 2.9.10): a CRAM and a BAM tumour/normal pair within a BED's regions,
and the CRAM pair with candidate indels and ``--exome`` ("best-practice").
As upstream, the calls are checked by nft-vcf's ``variantsMD5`` (the md5 of
the record lines) and the indices by name.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.strelka import StrelkaSomaticResult, strelka_somatic
from tests.framework import env, fixture

HS = "genomics/homo_sapiens/"
CHR21 = HS + "genome/chr21/sequence/"


async def strelka_variants_md5(vcf: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of a gzipped VCF's record lines joined without newlines, without leading zeros."""
    async with vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    records = [x for x in lines if not x.startswith("#")]
    return format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")


async def pair_inputs(alignments: str, index_suffix: str) -> dict[str, File]:
    normal, normal_index, tumor, tumor_index, fasta, fai, bed, bed_tbi = await asyncio.gather(
        fixture(alignments.format("test")),
        fixture(alignments.format("test") + index_suffix),
        fixture(alignments.format("test2")),
        fixture(alignments.format("test2") + index_suffix),
        fixture(CHR21 + "genome.fasta"),
        fixture(CHR21 + "genome.fasta.fai"),
        fixture(CHR21 + "multi_intervals.bed.gz"),
        fixture(CHR21 + "multi_intervals.bed.gz.tbi"),
    )
    return {
        "normal": normal,
        "normal_index": normal_index,
        "tumor": tumor,
        "tumor_index": tumor_index,
        "fasta": fasta,
        "fai": fai,
        "target_bed": bed,
        "target_bed_tbi": bed_tbi,
    }


async def check_calls(result: StrelkaSomaticResult, label: str, indels: str, snvs: str) -> None:
    assert result.vcf_indels.name == "test.somatic_indels.vcf.gz", f"{label}: {result.vcf_indels.name}"
    assert result.vcf_indels_tbi.name == "test.somatic_indels.vcf.gz.tbi", f"{label}: {result.vcf_indels_tbi.name}"
    assert result.vcf_snvs.name == "test.somatic_snvs.vcf.gz", f"{label}: {result.vcf_snvs.name}"
    assert result.vcf_snvs_tbi.name == "test.somatic_snvs.vcf.gz.tbi", f"{label}: {result.vcf_snvs_tbi.name}"
    actual = await strelka_variants_md5(result.vcf_indels)
    assert actual == indels, f"{label}: indels variantsMD5 {actual}, expected {indels}"
    actual = await strelka_variants_md5(result.vcf_snvs)
    assert actual == snvs, f"{label}: snvs variantsMD5 {actual}, expected {snvs}"


CRAM = HS + "illumina/cram/{}.paired_end.recalibrated.sorted.cram"
BAM = HS + "illumina/bam/{}.paired_end.recalibrated.sorted.bam"


@env.task
async def test_strelka_somatic_cram() -> None:
    # upstream case: strelka/somatic "human - cram"
    result = await strelka_somatic(**await pair_inputs(CRAM, ".crai"), prefix="test")
    await check_calls(result, "cram", "6b8a01b52186fcaf3ec74e5bf1b5426", "866236c526138e4e08fdc2fa87b08417")


@env.task
async def test_strelka_somatic_bam() -> None:
    # upstream case: strelka/somatic "human - bam"
    result = await strelka_somatic(**await pair_inputs(BAM, ".bai"), prefix="test")
    await check_calls(result, "bam", "6b8a01b52186fcaf3ec74e5bf1b5426", "866236c526138e4e08fdc2fa87b08417")


@env.task
async def test_strelka_somatic_cram_best_practice() -> None:
    # upstream case: strelka/somatic "human - cram - best-practice" (ext.args "--exome", with candidate indels)
    candidates = HS + "illumina/vcf/test.genome_21.somatic_sv.vcf.gz"
    inputs, indel_candidates, indel_candidates_tbi = await asyncio.gather(
        pair_inputs(CRAM, ".crai"), fixture(candidates), fixture(candidates + ".tbi")
    )
    result = await strelka_somatic(
        **inputs,
        prefix="test",
        indel_candidates=indel_candidates,
        indel_candidates_tbi=indel_candidates_tbi,
        args="--exome",
    )
    await check_calls(result, "best-practice", "a500521b4aa64e8f247d57fcb3195b39", "17d4e21e50940b4e75cca144af09ce66")


tests = [test_strelka_somatic_cram, test_strelka_somatic_bam, test_strelka_somatic_cram_best_practice]
