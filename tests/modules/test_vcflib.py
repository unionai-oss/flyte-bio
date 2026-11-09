"""Tests for flyte_bio.modules.vcflib.

The case is upstream's vcflib/vcffilter test (modules commit 9a1cc9f, vcflib
1.0.14), with its test config's ``-f "DP > 4"``. Upstream snapshots nft-vcf's
``variantsMD5`` (the md5 of the record lines), as this test does.
"""

import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.vcflib import vcflib_vcffilter
from tests.framework import env, fixture

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


tests = [test_vcflib_vcffilter_depth]
