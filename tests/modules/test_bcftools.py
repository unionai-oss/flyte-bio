"""Tests for flyte_bio.modules.bcftools.

Upstream's bcftools/stats snapshot records the report's first six lines,
which include the command line bcftools records.
"""

from flyte_bio.modules.bcftools import bcftools_stats
from tests.framework import env, fixture


@env.task
async def test_bcftools_stats() -> None:
    # upstream case: bcftools/stats "sarscov2 - vcf_gz"
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz")
    stats = await bcftools_stats(vcf, prefix="test")
    assert stats.path.endswith("/test.bcftools_stats.txt"), stats.path
    async with stats.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert lines[:6] == [
        "# This file was produced by bcftools stats (1.23.1+htslib-1.23.1) and can be plotted using plot-vcfstats.",
        "# The command line was:\tbcftools stats  test.vcf.gz",
        "#",
        "# Definition of sets:",
        "# ID\t[2]id\t[3]tab-separated file names",
        "ID\t0\ttest.vcf.gz",
    ], lines[:6]


tests = [test_bcftools_stats]
