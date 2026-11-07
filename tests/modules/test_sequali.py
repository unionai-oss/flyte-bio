"""Tests for flyte_bio.modules.sequali.

The cases are upstream's sequali tests (modules commit 9a1cc9f, Sequali
1.0.2). As upstream, the reports' names and the JSON report's ``summary``
are checked; the rest embeds the run's paths and date.
"""

import json

from flyte_bio.modules.sequali import SequaliResult, sequali
from tests.framework import env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


def summary(read_pair_info: str | None) -> dict:
    return {
        "mean_length": 138.97,
        "minimum_length": 72,
        "maximum_length": 151,
        "total_reads": 100,
        "q20_reads": 91,
        "total_bases": 13897,
        "q20_bases": 12922,
        "total_gc_bases": 5430,
        "total_n_bases": 0,
        "read_pair_info": read_pair_info,
    }


async def check(result: SequaliResult, expected: dict) -> None:
    assert (result.html.name, result.json.name) == ("test.html", "test.json"), (result.html.name, result.json.name)
    async with result.json.open("rb") as fh:
        report = json.loads(bytes(await fh.read()).decode())
    assert report["summary"] == expected, report["summary"]


@env.task
async def test_sequali_single_end() -> None:
    # upstream case: "sarscov2 - fastq single-end [sequali]"
    await check(await sequali(await fixture(FASTQ + "test_1.fastq.gz"), None, prefix="test"), summary(None))


@env.task
async def test_sequali_paired_end() -> None:
    # upstream case: "sarscov2 - fastq paired-end [sequali]"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    await check(await sequali(r1, r2, prefix="test"), summary("Read 1"))


@env.task
async def test_sequali_unaligned_bam() -> None:
    # upstream case: "sarscov2 - unaligned bam [sequali]"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.unaligned.bam")
    await check(await sequali(bam, None, prefix="test"), summary(None))


tests = [test_sequali_single_end, test_sequali_paired_end, test_sequali_unaligned_bam]
