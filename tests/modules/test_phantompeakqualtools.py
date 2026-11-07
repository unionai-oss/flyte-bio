"""Tests for flyte_bio.modules.phantompeakqualtools.

The cases and expected values are upstream's phantompeakqualtools tests
(modules commit 9a1cc9f, phantompeakqualtools 1.2.2): the metrics line's md5
and the plot's and R data's names. The metrics line records the BAM's name,
which the fixture keeps.
"""

import asyncio

from flyte_bio.modules.phantompeakqualtools import phantompeakqualtools
from tests.framework import assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/"


@env.task
async def test_phantompeakqualtools() -> None:
    # upstream cases: "sarscov2 - bam - single_end" and "sarscov2 - bam - paired_end"
    single, paired = await asyncio.gather(fixture(BAM + "test.single_end.bam"), fixture(BAM + "test.paired_end.bam"))
    results = await asyncio.gather(
        phantompeakqualtools(single, prefix="test"), phantompeakqualtools(paired, prefix="test")
    )
    for label, result, md5 in zip(
        ("single-end", "paired-end"),
        results,
        ("b01d976506b6fe45b66c821b1e8a1d15", "eed46e75eab119224f397a7a8b5924e6"),
        strict=True,
    ):
        assert (result.pdf.name, result.rdata.name) == ("test.spp.pdf", "test.spp.Rdata"), label
        await assert_md5(result.spp, md5, label=f"{label} test.spp.out")


tests = [test_phantompeakqualtools]
