"""Tests for flyte_bio.modules.vcftools (md5 from upstream's vcftools snapshot)."""

import asyncio

from flyte_bio.modules.vcftools import vcftools
from tests.framework import assert_md5, assert_nonempty, env, fixture


@env.task
async def test_vcftools() -> None:
    # upstream case: vcftools "sarscov2 - vcfgz" (--freq), plus the reports the pipeline asks for
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz")
    freq, count, qual, summary = await asyncio.gather(
        vcftools(vcf, prefix="test", args="--freq"),
        vcftools(vcf, prefix="test", args="--TsTv-by-count"),
        vcftools(vcf, prefix="test", args="--TsTv-by-qual"),
        vcftools(vcf, prefix="test", args="--FILTER-summary"),
    )
    frq = await freq.get_file("test.frq")
    assert frq is not None, "vcftools wrote no test.frq"
    await assert_md5(frq, "7f126655f17268fd1a338734f62868e9", label="test.frq")
    for results, name in ((count, "test.TsTv.count"), (qual, "test.TsTv.qual"), (summary, "test.FILTER.summary")):
        f = await results.get_file(name)
        assert f is not None, f"vcftools wrote no {name}"
        await assert_nonempty(f, label=name)


tests = [test_vcftools]
