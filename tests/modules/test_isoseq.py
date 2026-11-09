"""Tests for flyte_bio.modules.isoseq.

The case and expected md5s are upstream's isoseq/refine test (modules commit
9a1cc9f, isoseq 26.2.0), with its test config's prefix ``<id>.refine``.
nf-test hashes the BAM, its index and the CSV report as they are; the dataset
XML and the JSON summary are checked for content, as upstream. The BAM
header's ``@PG`` command line records ``-j``: upstream's snapshot is of a
4-thread run (reproduced with 4, not 1 or 2), so the test pins 4 CPUs.
"""

import asyncio

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.isoseq import isoseq_refine_cmd
from tests.framework import assert_md5, env, fixture

PACBIO = "genomics/homo_sapiens/pacbio/"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def output(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    assert found is not None, f"isoseq refine wrote no {name}"
    return found


@env.task
async def test_isoseq_refine() -> None:
    # upstream case: isoseq/refine "PacBio isoseq refine - Chimeras Removing test"
    bam, primers = await asyncio.gather(
        fixture(PACBIO + "bam/alz.ccs.fl.NEB_5p--NEB_Clontech_3p.bam"), fixture(PACBIO + "fasta/primers.fasta")
    )
    # The cache key ignores resources, so salt it to keep the 4-CPU run apart from another one.
    task = isoseq_refine_cmd.as_task().override(
        resources=flyte.Resources(cpu=4, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=4")
    )
    results = await task(bam=bam, primers=primers, prefix="test.refine", args="")
    expected = {
        "test.refine.bam": "e87acbb921049c0a506a5f7588f610be",
        "test.refine.bam.pbi": "5ed2cba7b7059040bbf0e3db48a39528",
        "test.refine.report.csv": "bbe3165430ec86f500cca980e22a8f82",
    }
    for name, md5 in expected.items():
        await assert_md5(await output(results, name), md5, label=name)
    xml = await text(await output(results, "test.refine.consensusreadset.xml"))
    assert "PacBio.DataSet.ConsensusReadSet" in xml, "no ConsensusReadSet in the XML"
    summary = await text(await output(results, "test.refine.filter_summary.report.json"))
    assert "pbcopper" in summary, "test.refine.filter_summary.report.json names no pbcopper"


tests = [test_isoseq_refine]
