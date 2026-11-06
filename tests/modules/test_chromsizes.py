"""Tests for flyte_bio.modules.chromsizes.

The md5s are the upstream module's snapshots from its last working version
(its current tests only check the deprecation failure). Sizes and indexes
are pure data, so they're compared byte for byte.
"""

from flyte_bio.modules.chromsizes import chrom_sizes
from tests.framework import assert_md5, env, fixture


@env.task
async def test_chrom_sizes() -> None:
    # upstream case: custom/getchromsizes "test_custom_getchromsizes"
    out = await chrom_sizes(await fixture("genomics/sarscov2/genome/genome.fasta"))
    assert out.sizes.path.endswith("/genome.fasta.sizes"), out.sizes.path
    assert out.fai.path.endswith("/genome.fasta.fai"), out.fai.path
    assert out.gzi is None
    await assert_md5(out.sizes, "a57c401f27ae5133823fb09fb21c8a3c", label="chromosome sizes")
    await assert_md5(out.fai, "9da2a56e2853dc8c0b86a9e7229c9fe5", label="fai")


@env.task
async def test_chrom_sizes_bgzip() -> None:
    # upstream case: custom/getchromsizes "test_custom_getchromsizes_bgzip"
    out = await chrom_sizes(await fixture("genomics/sarscov2/genome/genome.fasta.gz"))
    assert out.sizes.path.endswith("/genome.fasta.gz.sizes"), out.sizes.path
    assert out.gzi is not None and out.gzi.path.endswith("/genome.fasta.gz.gzi")
    await assert_md5(out.sizes, "a57c401f27ae5133823fb09fb21c8a3c", label="chromosome sizes")
    await assert_md5(out.fai, "9da2a56e2853dc8c0b86a9e7229c9fe5", label="fai")
    await assert_md5(out.gzi, "7dea362b3fac8e00956a4952a3d4f474", label="gzi")


tests = [test_chrom_sizes, test_chrom_sizes_bgzip]
