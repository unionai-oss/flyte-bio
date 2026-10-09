"""Tests for flyte_bio.modules.bedtools.

Each test is a `@env.task` that fetches fixtures, runs the wrapped
command, and asserts its output md5 matches the upstream snapshot value.
When upstream changes its snapshot, update the constant and re-run.
"""



from flyte_bio.modules.bedtools import bedtools_intersect, bedtools_sort
from tests.framework import assert_md5, env, fixture


@env.task
async def test_intersect_bed_bed() -> None:
    # upstream case: bedtools/intersect "sarscov2 - bed - bed"
    a = await fixture("genomics/sarscov2/genome/bed/test.bed")
    b = await fixture("genomics/sarscov2/genome/bed/test2.bed")
    out = await bedtools_intersect(a=a, b=[b])
    await assert_md5(
        out,
        expected="afcbf01c2f2013aad71dbe8e34f2c15c",
        label="bedtools intersect bed-bed",
    )


@env.task
async def test_sort_bed() -> None:
    # upstream case: bedtools/sort "test_bedtools_sort"
    i = await fixture("genomics/sarscov2/genome/bed/test.bed")
    out = await bedtools_sort(i=i)
    await assert_md5(
        out,
        expected="fe4053cf4de3aebbdfc3be2efb125a74",
        label="bedtools sort",
    )


tests = [test_intersect_bed_bed, test_sort_bed]
