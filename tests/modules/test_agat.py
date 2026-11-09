"""Tests for flyte_bio.modules.agat.

The cases and md5s are upstream's agat/convertgff2bed and agat/convertspgff2tsv
tests (modules commit 9a1cc9f, AGAT 1.6.1). Upstream's convertgff2bed test
config sets ``ext.args`` for a process named ``CONVERTGFF2BED``, which doesn't
match ``AGAT_CONVERTGFF2BED``, so its snapshot is of a run without args, as
here.
"""

from flyte_bio.modules.agat import agat_convertgff2bed, agat_convertspgff2tsv
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


tests = [test_agat_convertgff2bed, test_agat_convertspgff2tsv]
