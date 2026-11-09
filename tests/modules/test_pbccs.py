"""Tests for flyte_bio.modules.pbccs.

The case and expected md5s are upstream's pbccs test (modules commit 9a1cc9f,
ccs 6.4.0): chunk 2 of 3 of ``alz.bam``, with its test config's
``--min-rq 0.9``. nf-test hashes the BAM, its index and the text report as
they are; the JSON report and the metrics are checked for content, as
upstream. The BAM header's ``@PG`` command line records ``-j``: upstream's
snapshot is of a 4-thread run (reproduced with 4, not 1 or 2), so the test
pins 4 CPUs.
"""

import asyncio
import gzip

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.pbccs import pbccs_cmd
from tests.framework import assert_md5, env, fixture

PACBIO_BAM = "genomics/homo_sapiens/pacbio/bam/"


async def output(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    assert found is not None, f"ccs wrote no {name}"
    return found


@env.task
async def test_pbccs() -> None:
    # upstream case: pbccs "PacBio pbccs - Compute CSS test"
    bam, pbi = await asyncio.gather(fixture(PACBIO_BAM + "alz.bam"), fixture(PACBIO_BAM + "alz.bam.pbi"))
    # The cache key ignores resources, so salt it to keep the 4-CPU run apart from another one.
    task = pbccs_cmd.as_task().override(
        resources=flyte.Resources(cpu=4, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=4")
    )
    results = await task(bam=bam, pbi=pbi, chunk_num=2, chunk_on=3, prefix="alz", args="--min-rq 0.9")
    await assert_md5(
        await output(results, "alz.chunk2.bam"), "4f99cfdecba7d1d140530608afffa822", label="alz.chunk2.bam"
    )
    await assert_md5(
        await output(results, "alz.chunk2.bam.pbi"), "2a8630a2a6f58cb65ecb9c1c1ff705fb", label="alz.chunk2.bam.pbi"
    )
    await assert_md5(
        await output(results, "alz.chunk2.report.txt"),
        "67583febb144333e3f854c7bf4a17a50",
        label="alz.chunk2.report.txt",
    )
    async with (await output(results, "alz.chunk2.report.json")).open("rb") as fh:
        report = bytes(await fh.read()).decode()
    assert "pbcopper" in report, "alz.chunk2.report.json names no pbcopper"
    async with (await output(results, "alz.chunk2.metrics.json.gz")).open("rb") as fh:
        first_line = gzip.decompress(bytes(await fh.read())).decode().splitlines()[0]
    assert "zmws" in first_line, f"alz.chunk2.metrics.json.gz: first line has no zmws: {first_line[:200]}"


tests = [test_pbccs]
