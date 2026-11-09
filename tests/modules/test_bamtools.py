"""Tests for flyte_bio.modules.bamtools.

The cases and expected md5s are upstream's bamtools/stats, bamtools/split and
bamtools/convert tests (modules commit 9a1cc9f, BamTools 2.5.2); split runs
with its test config's ``-reference``, convert with each case's ``-format``.
Split's snapshot hashes the whole BAMs, which BamTools writes without a
command-line record, so they're compared byte for byte after unpacking the
archive. Upstream's convert error cases fail on a bad ``-format``; here the
wrapper raises before any task runs.
"""

import asyncio
import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bamtools import CONVERT_FORMAT_ERROR, bamtools_convert, bamtools_split, bamtools_stats
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


SARSCOV2_BAM = "genomics/sarscov2/illumina/bam/test.paired_end.bam"


async def check_convert_rejects(args: str) -> None:
    bam = await fixture(SARSCOV2_BAM)
    try:
        await bamtools_convert(bam, args=args, prefix="test")
    except ValueError as e:
        assert CONVERT_FORMAT_ERROR in str(e), str(e)
    else:
        raise AssertionError(f"bamtools convert accepted args {args!r}")


async def check_convert(fmt: str, md5: str) -> None:
    converted = await bamtools_convert(await fixture(SARSCOV2_BAM), args=f"-format {fmt}", prefix="test")
    assert converted.name == f"test.{fmt}", converted.name
    await assert_md5(converted, md5, label=f"test.{fmt}")


@env.task
async def test_bamtools_convert_ext_error() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_ext_error"
    await check_convert_rejects("-format vcf")


@env.task
async def test_bamtools_convert_noext_error() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_noext_error"
    await check_convert_rejects("")


@env.task
async def test_bamtools_convert_bed() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_bed"
    await check_convert("bed", "4e34cc15bf31e700f5f3a9f8fffb6c81")


@env.task
async def test_bamtools_convert_fasta() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_fasta"
    await check_convert("fasta", "52aeacf78571862b7e97c7d44ac8f827")


@env.task
async def test_bamtools_convert_fastq() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_fastq"
    await check_convert("fastq", "e591c48daad2c56638e5d6f21f1f71c5")


@env.task
async def test_bamtools_convert_json() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_json"
    await check_convert("json", "9c44279e6da864f30b52a79a14dcb2bd")


@env.task
async def test_bamtools_convert_pileup() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_pileup"
    await check_convert("pileup", "e5a3cb4a3e1bf980a575fafce6a2826f")


@env.task
async def test_bamtools_convert_sam() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_sam"
    await check_convert("sam", "6f559e53e1ea7ff3dc919c8c4695de96")


@env.task
async def test_bamtools_convert_yaml() -> None:
    # upstream case: bamtools/convert "test_bamtools_convert_yaml"
    await check_convert("yaml", "68b56f198da036fef33e150eb773dc3b")


tests = [
    test_bamtools_stats,
    test_bamtools_split_single_bam,
    test_bamtools_split_two_bams,
    test_bamtools_convert_ext_error,
    test_bamtools_convert_noext_error,
    test_bamtools_convert_bed,
    test_bamtools_convert_fasta,
    test_bamtools_convert_fastq,
    test_bamtools_convert_json,
    test_bamtools_convert_pileup,
    test_bamtools_convert_sam,
    test_bamtools_convert_yaml,
]
