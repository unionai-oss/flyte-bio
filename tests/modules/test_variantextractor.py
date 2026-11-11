"""Tests for flyte_bio.modules.variantextractor.

The case is upstream's variantextractor test (modules commit 9a1cc9f,
variant-extractor 5.1.0): the extracted VCF's md5.
"""

from flyte_bio.modules.variantextractor import variantextractor
from tests.framework import assert_md5, env, fixture


@env.task
async def test_variantextractor_vcf() -> None:
    vcf = await fixture("genomics/homo_sapiens/genome/vcf/ped/justhusky_minimal.vcf.gz")
    result = await variantextractor(vcf, prefix="test")
    assert result.vcf.name == "test.extracted.vcf", result.vcf.name
    await assert_md5(result.vcf, "150557d6bca69c57392ad17d363fc95c")


tests = [test_variantextractor_vcf]
