"""Tests for flyte_bio.modules.bedtools.

Each test is a `@env.task` that fetches fixtures, runs the wrapped
command, and asserts its output md5 matches the upstream snapshot value.
When upstream changes its snapshot, update the constant and re-run.
"""



from flyte_bio.modules.bedtools import (
    bedtools_bamtobed,
    bedtools_complement,
    bedtools_getfasta,
    bedtools_intersect,
    bedtools_maskfasta,
    bedtools_sort,
)
from tests.framework import assert_md5, env, fixture


@env.task
async def test_intersect_bed_bed() -> None:
    # upstream case: bedtools/intersect "sarscov2 - bed - bed"
    a = await fixture("genomics/sarscov2/genome/bed/test.bed")
    b = await fixture("genomics/sarscov2/genome/bed/test2.bed")
    out = await bedtools_intersect(a=a, b=[b])
    await assert_md5(
        out,
        expected="afcbf01c2f2013aad71dbe8e34f2c15c",
        label="bedtools intersect bed-bed",
    )


@env.task
async def test_sort_bed() -> None:
    # upstream case: bedtools/sort "test_bedtools_sort"
    i = await fixture("genomics/sarscov2/genome/bed/test.bed")
    out = await bedtools_sort(i=i)
    await assert_md5(
        out,
        expected="fe4053cf4de3aebbdfc3be2efb125a74",
        label="bedtools sort",
    )


@env.task
async def test_bedtools_getfasta() -> None:
    # upstream case: bedtools/getfasta "sarscov2 - bed - fasta"
    bed = await fixture("genomics/sarscov2/genome/bed/test.bed")
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    out = await bedtools_getfasta(bed, fasta, prefix="test")
    assert out.name == "test.fa", out.name
    await assert_md5(out, "41c3a45a57a16c04f828d8f8bb52df70", label="bedtools getfasta")


@env.task
async def test_bedtools_maskfasta() -> None:
    # upstream case: bedtools/maskfasta "test-bedtools-maskfasta"
    bed = await fixture("genomics/sarscov2/genome/bed/test.bed")
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    out = await bedtools_maskfasta(bed, fasta, prefix="test")
    assert out.name == "test.fa", out.name
    await assert_md5(out, "f4f6749698f11074228d2c79338e3b9c", label="bedtools maskfasta")


@env.task
async def test_bedtools_bamtobed_bam() -> None:
    # upstream case: bedtools/bamtobed "sarscov2 - bam" (modules commit 9a1cc9f)
    bam = await fixture("genomics/sarscov2/illumina/bam/test.single_end.bam")
    out = await bedtools_bamtobed(bam, prefix="test")
    assert out.name == "test.bed", out.name
    await assert_md5(out, "a6a299bd39dc56225f8029c05ea97dcb", label="bedtools bamtobed")


@env.task
async def test_bedtools_complement_bed() -> None:
    # upstream case: bedtools/complement "sarscov2" (modules commit 9a1cc9f)
    bed = await fixture("genomics/sarscov2/genome/bed/test.bed")
    sizes = await fixture("genomics/sarscov2/genome/genome.sizes")
    out = await bedtools_complement(bed, sizes, prefix="test_out")
    assert out.name == "test_out.bed", out.name
    await assert_md5(out, "d71df7e47aec0661c27b71e483e727f9", label="bedtools complement")


tests = [
    test_intersect_bed_bed,
    test_sort_bed,
    test_bedtools_getfasta,
    test_bedtools_maskfasta,
    test_bedtools_bamtobed_bam,
    test_bedtools_complement_bed,
]
