"""Tests for flyte_bio.modules.fastqdl.

The cases and expected md5s are upstream's fastqdl tests (modules commit
9a1cc9f, fastq-dl 4.0.1). They download about 320 MB from ENA, so they need
network. nf-test hashes the reads decompressed; the reads are too big for
the test pod to decompress in memory, so they're hashed streaming from disk.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.fastqdl import FastqdlResult, fastqdl
from tests.framework import assert_md5, env

READS = {
    "ERR5069949.fastq.gz": "c931a862c2126788432684cca99db416",
    "ERR5069949_1.fastq.gz": "b643cdf5f6a8fae5063fb2b7c172071d",
    "ERR5069949_2.fastq.gz": "685c1b66d8bf0c20a37402a0578acca2",
}
RUN_INFO_MD5 = "a24e465dd7c6e920afd443d743ebbb40"


def gunzipped_md5_of(path: str) -> str:
    digest = hashlib.md5()
    with gzip.open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


async def assert_streamed_gunzipped_md5(file: File, expected: str, *, label: str) -> None:
    """Like ``assert_gunzipped_md5``, without holding the file in memory."""
    actual = await asyncio.to_thread(gunzipped_md5_of, str(await file.download()))
    assert actual == expected, f"{label}: gunzipped md5 mismatch (expected {expected}, got {actual})"


async def check(result: FastqdlResult) -> None:
    names = [f.name for f in result.fastq]
    assert names == list(READS), f"fastq-dl wrote {names}, expected {list(READS)}"
    for read in result.fastq:
        await assert_streamed_gunzipped_md5(read, READS[read.name], label=read.name)
    assert result.runinfo.name == "test-run-info.tsv", f"run info named {result.runinfo.name!r}"
    await assert_md5(result.runinfo, RUN_INFO_MD5, label="test-run-info.tsv")
    assert result.runmergers is None, "fastq-dl merged runs without grouping"


@env.task
async def test_fastqdl_run() -> None:
    await check(await fastqdl("ERR5069949", prefix="test"))


@env.task
async def test_fastqdl_experiment() -> None:
    await check(await fastqdl("ERX4876079", prefix="test"))


tests = [test_fastqdl_run, test_fastqdl_experiment]
