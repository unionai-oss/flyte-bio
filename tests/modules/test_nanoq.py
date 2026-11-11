"""Tests for flyte_bio.modules.nanoq.

The cases and expected md5s are upstream's nanoq tests (modules commit
9a1cc9f, nanoq 0.10.0): nanopore reads written uncompressed and compressed
with gzip, bzip2 and lzma. nf-test hashes gzipped files decompressed, and
the others as they are.
"""

from flyte_bio.modules.nanoq import NanoqResult, nanoq
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

STATS_MD5 = "5ab32af3352dfeca8268e10edf6e4dbe"
READS_MD5 = "7567d853ada6ac142332619d0b541d76"


async def nanoq_test(output_format: str) -> NanoqResult:
    reads = await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz")
    result = await nanoq(reads, output_format, prefix="test_filtered")
    assert result.reads.name == f"test_filtered.{output_format}", result.reads.name
    await assert_md5(result.stats, STATS_MD5, label="test_filtered.stats")
    return result


@env.task
async def test_nanoq_uncompressed() -> None:
    # upstream case: nanoq "sarscov2 - nanopore_uncompressed"
    result = await nanoq_test("fastq")
    await assert_md5(result.reads, READS_MD5, label="test_filtered.fastq")


@env.task
async def test_nanoq_compressed_gz() -> None:
    # upstream case: nanoq "sarscov2 - nanopore_compressed_gz"
    result = await nanoq_test("fastq.gz")
    await assert_gunzipped_md5(result.reads, READS_MD5, label="test_filtered.fastq.gz")


@env.task
async def test_nanoq_compressed_bz2() -> None:
    # upstream case: nanoq "sarscov2 - nanopore_compressed_bz2"
    result = await nanoq_test("fastq.bz2")
    await assert_md5(result.reads, "b53cf14fd4eb5b16c459c41f03cc8a4b", label="test_filtered.fastq.bz2")


@env.task
async def test_nanoq_compressed_lzma() -> None:
    # upstream case: nanoq "sarscov2 - nanopore_compressed_lzma"
    result = await nanoq_test("fastq.lzma")
    await assert_md5(result.reads, "65dda701689f913734dc245b68c89e07", label="test_filtered.fastq.lzma")


tests = [test_nanoq_uncompressed, test_nanoq_compressed_gz, test_nanoq_compressed_bz2, test_nanoq_compressed_lzma]
