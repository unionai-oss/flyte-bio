"""Tests for flyte_bio.modules.bamtools.

The cases and expected md5s are upstream's bamtools/stats and bamtools/split
tests (modules commit 9a1cc9f, BamTools 2.5.2); split runs with its test
config's ``-reference``. Split's snapshot hashes the whole BAMs, which BamTools
writes without a command-line record, so they're compared byte for byte after
unpacking the archive.
"""

import asyncio
import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bamtools import bamtools_split, bamtools_stats
from tests.framework import assert_md5, env, fixture

HUMAN_BAM = "genomics/homo_sapiens/illumina/bam/"


def unpacked_md5s(archive: str) -> dict[str, str]:
    root = Path(tempfile.mkdtemp())
    with tarfile.open(archive) as tar:
        tar.extractall(root, filter="data")
    return {p.name: hashlib.md5(p.read_bytes()).hexdigest() for p in sorted(root.rglob("*.bam"))}


async def split_md5s(archive: File) -> dict[str, str]:
    return await asyncio.to_thread(unpacked_md5s, str(await archive.download()))


@env.task
async def test_bamtools_stats() -> None:
    # upstream case: bamtools/stats "test-bamtools-stats"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.bam")
    stats = await bamtools_stats(bam, prefix="test")
    assert stats.name == "test.bam.stats", stats.name
    await assert_md5(stats, "b7c3f01682673a652d664a2889b15b66", label="test.bam.stats")


@env.task
async def test_bamtools_split_single_bam() -> None:
    # upstream case: bamtools/split "single_bam"
    bam = await fixture(HUMAN_BAM + "test.paired_end.sorted.bam")
    archive = await bamtools_split([bam], args="-reference", prefix="test")
    assert archive.name == "test.split.tar.gz", archive.name
    expected = {
        "test.REF_chr22.bam": "b7dc50e0edf9c6bfc2e3b0e6d074dc07",
        "test.REF_unmapped.bam": "e0754bf72c51543b2d745d96537035fb",
    }
    actual = await split_md5s(archive)
    assert actual == expected, f"split BAMs: expected {expected}, got {actual}"


@env.task
async def test_bamtools_split_two_bams() -> None:
    # upstream case: bamtools/split "two_bams"
    bams = [
        await fixture(HUMAN_BAM + "test.paired_end.sorted.bam"),
        await fixture(HUMAN_BAM + "test2.paired_end.sorted.bam"),
    ]
    archive = await bamtools_split(bams, args="-reference", prefix="test")
    expected = {
        "test.REF_chr22.bam": "585675bea34c48ebe9db06a561d4b4fa",
        "test.REF_unmapped.bam": "16ad644c87b9471f3026bc87c98b4963",
    }
    actual = await split_md5s(archive)
    assert actual == expected, f"split BAMs: expected {expected}, got {actual}"


tests = [test_bamtools_stats, test_bamtools_split_single_bam, test_bamtools_split_two_bams]
