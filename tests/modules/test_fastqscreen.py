"""Tests for flyte_bio.modules.fastqscreen.

The cases and expected md5s are upstream's fastqscreen/buildfromindex and
fastqscreen/fastqscreen tests (modules commit 9a1cc9f, FastQ Screen 0.16.0).
As upstream, the database is built from Bowtie 2 indexes of the SARS-CoV-2
genome and human chr21; the "tags" case runs with upstream's test config
(``--tag --filter 10``). nf-test computes the tagged reads' md5s
decompressed. The module's 4 CPUs match upstream's test CPUs.
"""

import asyncio

from flyte.io import Dir, File

from flyte_bio.modules.bowtie2 import bowtie2_build
from flyte_bio.modules.fastqscreen import fastqscreen, fastqscreen_database
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

SARSCOV2 = "genomics/sarscov2/genome/genome.fasta"
HUMAN = "genomics/homo_sapiens/genome/chr21/sequence/genome.fasta"
READS = "genomics/sarscov2/illumina/fastq/test_1.fastq.gz"

INDEXES = {
    "dir1": {
        "genome.1.bt2": "cbe3d0bbea55bc57c99b4bfa25b5fbdf",
        "genome.2.bt2": "47b153cd1319abc88dda532462651fcf",
        "genome.3.bt2": "4ed93abba181d8dfab2e303e33114777",
        "genome.4.bt2": "c25be5f8b0378abf7a58c8a880b87626",
        "genome.rev.1.bt2": "52be6950579598a990570fbcf5372184",
        "genome.rev.2.bt2": "e3b4ef343dea4dd571642010a7d09597",
    },
    "dir2": {
        "genome.1.bt2": "2fbc8eeaf480f03b3a2362a782fa5755",
        "genome.2.bt2": "a5300e3cf590d4cd8bd5521cac6337d6",
        "genome.3.bt2": "d2a7a0e1c9891e847a15781c0695b84e",
        "genome.4.bt2": "a55dc4b09e3e00586de6963800c95b9a",
        "genome.rev.1.bt2": "93c07f7d04e5da6bf1322ad5f46ac806",
        "genome.rev.2.bt2": "c3fa30132c0c9b9bc7868059e8e23466",
    },
}

HTML_HEAD = [
    "<!DOCTYPE html>",
    "<html>",
    "",
    "",
    "<!-- Style, and Plotly code -->",
    "<head>",
    '\t<meta http-equiv="content-type" content="text/html; charset=UTF-8">',
    "\t<title>FastQ Screen Processing Report - test_1.fastq.gz</title>",
    '\t<style type="text/css">',
    "\t\tbody {",
    "\t\t\tfont-family: Arial, sans-serif;",
]


async def database() -> Dir:
    fastas = await asyncio.gather(fixture(SARSCOV2), fixture(HUMAN))
    indexes = await asyncio.gather(*(bowtie2_build(fasta) for fasta in fastas))
    return await fastqscreen_database(["sarscov2", "human"], list(indexes))


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_fastqscreen_buildfromindex() -> None:
    # upstream case: fastqscreen/buildfromindex "sarscov2 - human"
    db = await database()
    conf = await db.get_file("fastq_screen.conf")
    assert conf is not None, "database has no fastq_screen.conf"
    await assert_md5(conf, "ca56e866c7006a46bdfe4e751e787265", label="fastq_screen.conf")
    for folder, md5s in INDEXES.items():
        for name, md5 in md5s.items():
            file = await db.get_file(f"{folder}/{name}")
            assert file is not None, f"database has no {folder}/{name}"
            await assert_md5(file, md5, label=f"{folder}/{name}")


@env.task
async def test_fastqscreen_screen() -> None:
    # upstream case: fastqscreen/fastqscreen "sarscov2 - human"
    result = await fastqscreen([await fixture(READS)], await database())
    names = [f.name for f in result.txt + result.png + result.html]
    assert names == ["test_1_screen.txt", "test_1_screen.png", "test_1_screen.html"], names
    assert result.fastq == [], [f.name for f in result.fastq]


@env.task
async def test_fastqscreen_tags() -> None:
    # upstream case: fastqscreen/fastqscreen "sarscov2 - human - tags"
    result = await fastqscreen([await fixture(READS)], await database(), args="--tag --filter 10")
    assert [f.name for f in result.txt] == ["test_1_screen.txt"]
    await assert_md5(result.txt[0], "7c0e31428663cce39b49f50cff70b02d", label="test_1_screen.txt")
    names = [f.name for f in result.fastq]
    assert names == ["test_1.tagged.fastq.gz", "test_1.tagged_filter.fastq.gz"], names
    await assert_gunzipped_md5(result.fastq[0], "f742b162c43ce28f80b89608d5c47f3d", label=names[0])
    await assert_gunzipped_md5(result.fastq[1], "28527a76bb0bb3fce0ee76afe01e90aa", label=names[1])
    head = (await text(result.html[0])).splitlines()[:11]  # as Groovy's readLines, without the \r
    assert head == HTML_HEAD, head
    assert [f.name for f in result.png] == ["test_1_screen.png"]


tests = [test_fastqscreen_buildfromindex, test_fastqscreen_screen, test_fastqscreen_tags]
