"""Tests for flyte_bio.modules.mosdepth (md5s from upstream's mosdepth snapshot).

nf-test hashes ``.gz`` outputs after decompressing them.
"""

import asyncio

from flyte_bio.modules.mosdepth import mosdepth
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"


@env.task
async def test_mosdepth_bam() -> None:
    # upstream case: mosdepth "homo_sapiens - bam, bai, []"
    bam, bai = await asyncio.gather(fixture(BAM), fixture(BAM + ".bai"))
    results = await mosdepth(bam, bai, prefix="test")
    for name, md5 in {
        "test.mosdepth.global.dist.txt": "e82e90c7d508a135b5a8a7cd6933452e",
        "test.mosdepth.summary.txt": "4f0d231060cbde4efdd673863bd2fb59",
    }.items():
        f = await results.get_file(name)
        assert f is not None, f"mosdepth wrote no {name}"
        await assert_md5(f, md5, label=name)
    per_base = await results.get_file("test.per-base.bed.gz")
    assert per_base is not None
    await assert_gunzipped_md5(per_base, "da6db0fb375a3053a89db8c935eebbaa", label="per-base.bed.gz")


@env.task
async def test_mosdepth_bed() -> None:
    # upstream case: mosdepth "homo_sapiens - bam, bai, bed"
    bam, bai, bed = await asyncio.gather(
        fixture(BAM), fixture(BAM + ".bai"), fixture("genomics/homo_sapiens/genome/genome.bed")
    )
    results = await mosdepth(bam, bai, prefix="test", bed=bed)
    for name, md5 in {
        "test.mosdepth.global.dist.txt": "e82e90c7d508a135b5a8a7cd6933452e",
        "test.mosdepth.summary.txt": "96c037f769974b904beb53edc4f56d82",
        "test.mosdepth.region.dist.txt": "e82e90c7d508a135b5a8a7cd6933452e",
    }.items():
        f = await results.get_file(name)
        assert f is not None, f"mosdepth wrote no {name}"
        await assert_md5(f, md5, label=name)


tests = [test_mosdepth_bam, test_mosdepth_bed]
