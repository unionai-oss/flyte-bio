"""Tests for flyte_bio.modules.vcflib.

The case is upstream's vcflib/vcffilter test (modules commit 9a1cc9f, vcflib
1.0.14), with its test config's ``-f "DP > 4"``. Upstream snapshots nft-vcf's
``variantsMD5`` (the md5 of the record lines), as this test does.

The vcffixup cases are upstream's vcflib/vcffixup tests at the same commit:
the outputs' names, and the bgzipped VCF's ``variantsMD5``.

The vcfuniq case is upstream's vcflib/vcfuniq test at the same commit, which
snapshots the output's decompressed md5 (bgzip writes no dates or names).
"""

import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.vcflib import vcflib_vcffilter, vcflib_vcffixup, vcflib_vcfuniq
from tests.framework import assert_gunzipped_md5, env, fixture

VCF = "genomics/sarscov2/illumina/vcf/test.vcf.gz"


async def vcflib_variants_md5(file: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of a gzipped VCF's record lines, joined without newlines."""
    async with file.open("rb") as fh:
        text = gzip.decompress(bytes(await fh.read())).decode()
    body = "".join(line for line in text.splitlines() if not line.startswith("#"))
    return format(int(hashlib.md5(body.encode()).hexdigest(), 16), "x")


@env.task
async def test_vcflib_vcffilter_depth() -> None:
    # upstream case: vcflib/vcffilter "test-vcflib-vcffilter"
    result = await vcflib_vcffilter(
        await fixture(VCF), info_filters=["DP > 4"], tbi=await fixture(VCF + ".tbi"), prefix="test.filter"
    )
    assert result.name == "test.filter.vcf.gz", f"wrote {result.name}"
    variants = await vcflib_variants_md5(result)
    assert variants == "86a26a6b691c05bc7f6e74772b44543d", f"variantsMD5 {variants}"


@env.task
async def test_vcflib_vcffixup() -> None:
    # upstream cases: vcflib/vcffixup "test-vcflib-vcffixup - vcf.gz - no index" and "... - bcf.gz - csi"
    vcf, tbi = await fixture(VCF), await fixture(VCF + ".tbi")
    as_vcf = await vcflib_vcffixup(vcf, tbi=tbi, prefix="test.fixed", args2="-Oz")
    assert as_vcf.vcf.name == "test.fixed.vcf.gz" and as_vcf.index is None, f"wrote {as_vcf.vcf.name}"
    variants = await vcflib_variants_md5(as_vcf.vcf)
    assert variants == "129354635f189827c91335aab65b084", f"variantsMD5 {variants}"
    as_bcf = await vcflib_vcffixup(vcf, tbi=tbi, prefix="test.fixed", args2="-Ob --write-index=csi")
    names = (as_bcf.vcf.name, as_bcf.index.name if as_bcf.index else None)
    assert names == ("test.fixed.bcf.gz", "test.fixed.bcf.gz.csi"), f"wrote {names}"


@env.task
async def test_vcflib_vcfuniq() -> None:
    # upstream case: vcflib/vcfuniq "test-vcflib-vcfuniq"
    result = await vcflib_vcfuniq(await fixture(VCF), tbi=await fixture(VCF + ".tbi"), prefix="test.uniq")
    assert result.name == "test.uniq.vcf.gz", f"wrote {result.name}"
    await assert_gunzipped_md5(result, "8e722884ffb75155212a3fc053918766", label="test.uniq.vcf.gz")


tests = [test_vcflib_vcffilter_depth, test_vcflib_vcffixup, test_vcflib_vcfuniq]
