"""Tests for flyte_bio.modules.deepvariant.

The cases and expected md5s are upstream's deepvariant/rundeepvariant tests
(modules commit 9a1cc9f, DeepVariant 1.10.0), all with the WGS model and the
prefix ``test_out``: the homo_sapiens test BAM within ``chr22:0-40001``
against the plain and the bgzipped reference, the test CRAM within
``genome.bed``, and the CRAM again with chr22 haploid and
``genome.blacklist_intervals.bed`` as its pseudo-autosomal regions.

nf-test hashes the ``.vcf.gz`` outputs decompressed and the ``.tbi`` indexes
as they are. The VCFs, gVCFs and VCF indexes match upstream's md5s; the gVCF
indexes don't (see :func:`check`), so they get a content check.

Upstream quotes its region (``--regions="chr22:0-40001"``), which its shell
strips; ``args`` here isn't parsed by a shell, so it has no quotes.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.deepvariant import DeepVariantResult, deepvariant_rundeepvariant
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"
CRAM = "genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram"

DIPLOID = {
    "vcf": "d6ec6971b53fc10bf9ed23ddf6f544c1",
    "vcf_index": "efb861c33cf51739a4ae13a2e81126b6",
    "gvcf": "71bf2757abd13abed78479887f9de8f6",
    "gvcf_index": "90732f1ea2188c12fc4dab62da039f63",  # not reproduced, see check()
}
HAPLOID_CHR22 = {
    "vcf": "44664d7a287819aa3c18d5bd1aafed0b",
    "vcf_index": "7f57f737784a7102360ae06a636eac14",
    "gvcf": "61246dde2499b2301429eff9e96714dd",
    "gvcf_index": "3cb9cc54bbcd6bf44d2b835aeeaa5f9e",  # not reproduced, see check()
}


async def check(result: DeepVariantResult, expected: dict[str, str]) -> None:
    files: dict[str, File] = {
        "vcf": result.vcf,
        "vcf_index": result.vcf_index,
        "gvcf": result.gvcf,
        "gvcf_index": result.gvcf_index,
    }
    names = {key: file.name for key, file in files.items()}
    assert names == {
        "vcf": "test_out.vcf.gz",
        "vcf_index": "test_out.vcf.gz.tbi",
        "gvcf": "test_out.g.vcf.gz",
        "gvcf_index": "test_out.g.vcf.gz.tbi",
    }, names
    for key in ("vcf", "vcf_index", "gvcf"):
        assert_file_md5 = assert_gunzipped_md5 if key in ("vcf", "gvcf") else assert_md5
        await assert_file_md5(files[key], expected[key], label=files[key].name)
    # The gVCF's index records the gVCF's compressed block offsets, and our
    # gVCF compresses to different bytes than upstream's (its content matches):
    # check it's a tabix index instead of its md5.
    async with result.gvcf_index.open("rb") as fh:
        index = gzip.decompress(bytes(await fh.read()))
    assert index[:4] == b"TBI\x01", f"{result.gvcf_index.name} isn't a tabix index"
    assert result.report is None, result.report


@env.task
async def test_deepvariant_rundeepvariant_bam() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(BAM),
        await fixture(BAM + ".bai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        prefix="test_out",
        args="--regions=chr22:0-40001 --model_type=WGS",
    )
    await check(result, DIPLOID)


@env.task
async def test_deepvariant_rundeepvariant_cram_intervals() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(CRAM),
        await fixture(CRAM + ".crai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        intervals=await fixture(GENOME + "genome.bed"),
        prefix="test_out",
        args="--model_type=WGS",
    )
    await check(result, DIPLOID)


@env.task
async def test_deepvariant_rundeepvariant_cram_intervals_par_bed() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(CRAM),
        await fixture(CRAM + ".crai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        intervals=await fixture(GENOME + "genome.bed"),
        par_bed=await fixture(GENOME + "genome.blacklist_intervals.bed"),
        prefix="test_out",
        args="--model_type=WGS --haploid_contigs chr22",
    )
    await check(result, HAPLOID_CHR22)


@env.task
async def test_deepvariant_rundeepvariant_bam_fasta_gz() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(BAM),
        await fixture(BAM + ".bai"),
        await fixture(GENOME + "genome.fasta.gz"),
        await fixture(GENOME + "genome.fasta.gz.fai"),
        gzi=await fixture(GENOME + "genome.fasta.gz.gzi"),
        prefix="test_out",
        args="--regions=chr22:0-40001 --model_type=WGS",
    )
    await check(result, DIPLOID)


tests = [
    test_deepvariant_rundeepvariant_bam,
    test_deepvariant_rundeepvariant_cram_intervals,
    test_deepvariant_rundeepvariant_cram_intervals_par_bed,
    test_deepvariant_rundeepvariant_bam_fasta_gz,
]
