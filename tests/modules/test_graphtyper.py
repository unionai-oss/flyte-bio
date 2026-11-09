"""Tests for flyte_bio.modules.graphtyper.

The cases are upstream's graphtyper/genotype (GraphTyper 2.7.2) and
graphtyper/vcfconcatenate (GraphTyper 2.7.7) tests at modules commit 9a1cc9f.
VCF headers record the date and command line, so upstream snapshots nft-vcf's
``variantsMD5`` (the md5 of the record lines; the region file's case hashes
only the first VCF, as the second is empty), the indexes' md5s or names, and,
for the two-BAM case, only the output names; these tests assert the same.
"""

import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.graphtyper import graphtyper_genotype, graphtyper_vcfconcatenate
from tests.framework import assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/"
GENOME = "genomics/sarscov2/genome/"


async def graphtyper_variants_md5(file: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of a gzipped VCF's record lines, joined without newlines."""
    async with file.open("rb") as fh:
        text = gzip.decompress(bytes(await fh.read())).decode()
    body = "".join(line for line in text.splitlines() if not line.startswith("#"))
    return format(int(hashlib.md5(body.encode()).hexdigest(), 16), "x")


async def reference() -> tuple[File, File]:
    return await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.fasta.fai")


@env.task
async def test_graphtyper_genotype_single_region_file() -> None:
    # upstream case: graphtyper/genotype "single - region_file"
    ref, fai = await reference()
    result = await graphtyper_genotype(
        await fixture(BAM + "test.paired_end.sorted.bam"),
        await fixture(BAM + "test.paired_end.sorted.bam.bai"),
        ref,
        fai,
        region_file=await fixture(GENOME + "graphtyper/regions.txt"),
    )
    names = [f.name for f in result.vcf]
    assert names == ["000000001-000029829.vcf.gz", "000001000-000005000.vcf.gz"], f"vcfs: {names}"
    variants = await graphtyper_variants_md5(result.vcf[0])
    assert variants == "e238475710bfaafe703704b006568720", f"variantsMD5 {variants}"
    tbi = {
        "000000001-000029829.vcf.gz.tbi": "6cf39cb5edec1dad83e9bd85b68ae9d9",
        "000001000-000005000.vcf.gz.tbi": "4cb176febbc8c26d717a6c6e67b9c905",
    }
    assert [f.name for f in result.tbi] == list(tbi), f"indexes: {[f.name for f in result.tbi]}"
    for file in result.tbi:
        await assert_md5(file, tbi[file.name], label=file.name)


@env.task
async def test_graphtyper_genotype_multiple_region_args() -> None:
    # upstream case: graphtyper/genotype "multiple - region_args"
    ref, fai = await reference()
    result = await graphtyper_genotype(
        [await fixture(BAM + "test.paired_end.sorted.bam"), await fixture(BAM + "test.single_end.sorted.bam")],
        [await fixture(BAM + "test.paired_end.sorted.bam.bai"), await fixture(BAM + "test.single_end.sorted.bam.bai")],
        ref,
        fai,
        args="--region MT192765.1",
    )
    assert [f.name for f in result.vcf] == ["000000001-000029829.vcf.gz"], f"vcfs: {[f.name for f in result.vcf]}"
    assert [f.name for f in result.tbi] == ["000000001-000029829.vcf.gz.tbi"], f"tbi: {[f.name for f in result.tbi]}"


@env.task
async def test_graphtyper_vcfconcatenate_two_vcfs() -> None:
    # upstream case: graphtyper/vcfconcatenate "test-graphtyper-vcfconcatenate"
    vcfs = [
        await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz"),
        await fixture("genomics/sarscov2/illumina/vcf/test2.vcf.gz"),
    ]
    result = await graphtyper_vcfconcatenate(vcfs, prefix="test_all")
    assert result.vcf.name == "test_all.vcf.gz", f"wrote {result.vcf.name}"
    variants = await graphtyper_variants_md5(result.vcf)
    assert variants == "93f741cc7a648675dfc13f0a42a39788", f"variantsMD5 {variants}"
    assert result.tbi.name == "test_all.vcf.gz.tbi", f"wrote {result.tbi.name}"


tests = [
    test_graphtyper_genotype_single_region_file,
    test_graphtyper_genotype_multiple_region_args,
    test_graphtyper_vcfconcatenate_two_vcfs,
]
