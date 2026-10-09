"""Tests for flyte_bio.modules.htslib (md5s from upstream's htslib/bgziptabix snapshot).

nf-test hashes ``.gz`` outputs after decompressing them; the ``.tbi`` raw.
"""

from flyte_bio.modules.htslib import htslib_bgziptabix
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture


@env.task
async def test_bgziptabix_vcf() -> None:
    # upstream case: htslib/bgziptabix "sarscov2 - vcf - compress - index"
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf")
    gz, tbi = await htslib_bgziptabix(vcf, prefix="test", out_ext="vcf")
    assert gz.path.endswith("/test.vcf.gz") and tbi is not None, (gz.path, tbi)
    await assert_gunzipped_md5(gz, "8e722884ffb75155212a3fc053918766", label="test.vcf.gz")
    await assert_md5(tbi, "7f005943c935f2b55ba3f9d4802aa09f", label="test.vcf.gz.tbi")


@env.task
async def test_bgzip_no_index() -> None:
    # upstream case: htslib/bgziptabix "sarscov2 - gzip - (re)compress - no index" (a gzipped input)
    fastq = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    gz, tbi = await htslib_bgziptabix(fastq, prefix="test", out_ext="fastq", make_index=False)
    assert tbi is None
    await assert_gunzipped_md5(gz, "4161df271f9bfcd25d5845a1e220dbec", label="test.fastq.gz")


tests = [test_bgziptabix_vcf, test_bgzip_no_index]
