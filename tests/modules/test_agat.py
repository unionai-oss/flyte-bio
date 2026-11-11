"""Tests for flyte_bio.modules.agat.

The cases and md5s are upstream's agat/convertgff2bed, agat/convertspgff2tsv,
agat/spkeeplongestisoform, agat/convertspgxf2gxf and agat/spstatistics tests
(modules commit 9a1cc9f, AGAT 1.6.1). Upstream's convertgff2bed test
config sets ``ext.args`` for a process named ``CONVERTGFF2BED``, which doesn't
match ``AGAT_CONVERTGFF2BED``, so its snapshot is of a run without args, as
here.
"""

from flyte_bio.modules.agat import (
    agat_convertgff2bed,
    agat_convertspgff2tsv,
    agat_convertspgxf2gxf,
    agat_spkeeplongestisoform,
    agat_spstatistics,
)
from tests.framework import assert_md5, env, fixture

GFF3 = "genomics/sarscov2/genome/genome.gff3"


@env.task
async def test_agat_convertgff2bed() -> None:
    # upstream case: agat/convertgff2bed "sarscov2 genome - gff3"
    bed = await agat_convertgff2bed(await fixture(GFF3), prefix="test")
    assert bed.name == "test.bed", bed.name
    await assert_md5(bed, "d424ee1b9d9982f129412f491fca709e", label="bed")


@env.task
async def test_agat_convertspgff2tsv() -> None:
    # upstream case: agat/convertspgff2tsv "sarscov2 - genome [gff3]"
    tsv = await agat_convertspgff2tsv(await fixture(GFF3), prefix="test")
    assert tsv.name == "test.tsv", tsv.name
    await assert_md5(tsv, "3bd9dd9d17073d014bb97626cea9485a", label="tsv")


@env.task
async def test_agat_spkeeplongestisoform_gff() -> None:
    # upstream case: agat/spkeeplongestisoform "homo_sapiens - gff"
    gff = await fixture("genomics/eukaryotes/actinidia_chinensis/genome/chr1/genome.gff3.gz")
    longest = await agat_spkeeplongestisoform(gff, prefix="test")
    assert longest.name == "test.longest.gff", longest.name
    await assert_md5(longest, "019d7435f6c64eda77d290ce6393f3bc", label="longest.gff")


@env.task
async def test_agat_spkeeplongestisoform_gtf() -> None:
    # upstream case: agat/spkeeplongestisoform "homo_sapiens - gtf"
    longest = await agat_spkeeplongestisoform(await fixture("genomics/homo_sapiens/genome/genome.gtf"), prefix="test")
    await assert_md5(longest, "2a2a92171bb43592c1928393e6387540", label="longest.gff")


@env.task
async def test_agat_convertspgxf2gxf() -> None:
    # upstream case: agat/convertspgxf2gxf "sarscov2 genome [gtf]"
    gff = await agat_convertspgxf2gxf(await fixture("genomics/sarscov2/genome/genome.gtf"), prefix="test")
    assert gff.name == "test.agat.gff", gff.name
    await assert_md5(gff, "7d7e9bcd82a2f0bb7d8a38f85e82f0bc", label="agat.gff")


@env.task
async def test_agat_spstatistics() -> None:
    # upstream case: agat/spstatistics "sarscov2 genome [gff]"
    stats = await agat_spstatistics(await fixture(GFF3), prefix="test")
    assert stats.name == "test.stats.txt", stats.name
    await assert_md5(stats, "49b9c37029287f60be88a209abd2f58d", label="stats.txt")


tests = [
    test_agat_convertgff2bed,
    test_agat_convertspgff2tsv,
    test_agat_spkeeplongestisoform_gff,
    test_agat_spkeeplongestisoform_gtf,
    test_agat_convertspgxf2gxf,
    test_agat_spstatistics,
]
