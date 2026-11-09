"""Tests for flyte_bio.modules.pbcpgtools.

The cases and expected md5s are upstream's pbcpgtools/alignedbamtocpgscores
tests (modules commit 9a1cc9f, pb-CpG-tools 3.0.0): a downsampled PacBio BAM,
not phased (combined scores only) and phased (combined and per-haplotype
scores). nf-test hashes the ``.bed.gz`` files decompressed and the indexes
and bigWigs as they are.

The BED headers record the command line, ``--threads`` included, and
upstream's snapshot was made with 2 threads, so the tests pin 2 CPUs,
salting the cache (which ignores resources).
"""

import flyte
from flyte.io import File

from flyte_bio.modules.pbcpgtools import (
    CpgScores,
    PbcpgtoolsResult,
    pbcpgtools_alignedbamtocpgscores_cmd,
    pbcpgtools_scores,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

BAM = "genomics/homo_sapiens/pacbio/bam/NA037562_downsampled.pbmm2.repeats"
# Upstream's snapshot threads; the salt keeps these runs' cache entries apart from default-resource ones.
TWO_CPUS = {"resources": flyte.Resources(cpu=2, memory="8Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=2")}


async def pbcpgtools_two_cpus(bam: File, bai: File) -> PbcpgtoolsResult:
    task = pbcpgtools_alignedbamtocpgscores_cmd.as_task().override(**TWO_CPUS)
    return await pbcpgtools_scores(await task(bam=bam, bai=bai, prefix="test", args=""), "test")


async def check(scores: CpgScores | None, track: str, md5s: tuple[str, str, str]) -> None:
    assert scores is not None, f"no {track} scores"
    files = (scores.bed, scores.bed_index, scores.bigwig)
    names = [f.name for f in files]
    assert names == [f"test.{track}.bed.gz", f"test.{track}.bed.gz.tbi", f"test.{track}.bw"], names
    await assert_gunzipped_md5(scores.bed, md5s[0], label=scores.bed.name)
    await assert_md5(scores.bed_index, md5s[1], label=scores.bed_index.name)
    await assert_md5(scores.bigwig, md5s[2], label=scores.bigwig.name)


@env.task
async def test_pbcpgtools_alignedbamtocpgscores_notphased() -> None:
    result = await pbcpgtools_two_cpus(await fixture(BAM + ".bam"), await fixture(BAM + ".bam.bai"))
    await check(
        result.combined,
        "combined",
        ("bacdb71f8eec683930c6103f020277cd", "9ca5efde7981fec4e6a79b0b08999186", "77a572a00645df6743459805bd397cb9"),
    )
    assert result.hap1 is None and result.hap2 is None, "unphased reads gave haplotype scores"


@env.task
async def test_pbcpgtools_alignedbamtocpgscores_phased() -> None:
    result = await pbcpgtools_two_cpus(await fixture(BAM + ".phased.bam"), await fixture(BAM + ".phased.bam.bai"))
    await check(
        result.combined,
        "combined",
        ("861084b58d048a2a957a86068ec39bcf", "4325f78aeb593d00de8b56cffbc77712", "77a572a00645df6743459805bd397cb9"),
    )
    await check(
        result.hap1,
        "hap1",
        ("db37b440b02203f38c0d7c20b3dd62f3", "b9813e40b7d3365fb312ad35ac117bab", "df389b5415312c1da05a762177f7b50e"),
    )
    await check(
        result.hap2,
        "hap2",
        ("b5b7a4539d3906bd938010644b779196", "0213f587a323929cd2180862b01e8044", "b87e364ef8cf80790ddef4963fa5e55c"),
    )


tests = [
    test_pbcpgtools_alignedbamtocpgscores_notphased,
    test_pbcpgtools_alignedbamtocpgscores_phased,
]
