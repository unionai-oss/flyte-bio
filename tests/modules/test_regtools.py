"""Tests for flyte_bio.modules.regtools.

The cases and expected md5s are upstream's regtools/junctionsextract tests
(modules commit 9a1cc9f, regtools 1.0.0).
"""

from flyte.io import File

from flyte_bio.modules.regtools import regtools_junctionsextract
from tests.framework import assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.chr6.bam"


async def junctions(prefix: str, strand_specificity: str) -> File:
    return await regtools_junctionsextract(
        await fixture(BAM), await fixture(BAM + ".bai"), prefix=prefix, strand_specificity=strand_specificity
    )


@env.task
async def test_regtools_junctionsextract() -> None:
    # upstream case: regtools/junctionsextract "homo_sapiens - bam-bai" (an empty strandedness: XS)
    await assert_md5(await junctions("test", ""), "abcb59fd5c52fb16acf6ecc6528d26b7", label="junc")


@env.task
async def test_regtools_junctionsextract_first_strand() -> None:
    # upstream case: regtools/junctionsextract "homo_sapiens - bam-bai - first-strand"
    await assert_md5(await junctions("test", "RF"), "c0ccba2d76e9cb1588ba6d2e18c0742c", label="junc")


tests = [test_regtools_junctionsextract, test_regtools_junctionsextract_first_strand]
