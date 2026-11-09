"""Tests for flyte_bio.modules.humid.

The cases and expected md5s are upstream's humid tests (modules commit
9a1cc9f, HUMID 1.0.4): test reads trimmed by fastp with their UMIs moved
into the read names (``--umi --umi_loc=per_read --umi_len=8``), as in
upstream's setup, then deduplicated with statistics and annotated reads
(``-s -a``). nf-test hashes the gzipped reads decompressed, so these tests
do too. The log isn't snapshotted upstream, which checks that it reports
reading the data, as these tests do.
"""

from flyte.io import Dir, File

from flyte_bio.modules.fastp import fastp
from flyte_bio.modules.humid import HumidResult, humid
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
FASTP_ARGS = "--umi --umi_loc=per_read --umi_len=8"
HUMID_ARGS = "-s -a"


async def check_gunzipped(files: list[File], expected: dict[str, str], label: str) -> None:
    assert [f.name for f in files] == list(expected), f"{label}: {[f.name for f in files]}"
    for file in files:
        await assert_gunzipped_md5(file, expected[file.name], label=file.name)


async def check_stats(stats: Dir | None, expected: dict[str, str]) -> None:
    assert stats is not None, "humid wrote no statistics"
    files = {f.name: f for f in await stats.list_files()}
    assert sorted(files) == sorted(expected), f"statistics: {sorted(files)}"
    for name, md5 in expected.items():
        await assert_md5(files[name], md5, label=name)


async def check_log(result: HumidResult) -> None:
    assert result.log.name == "test.log", result.log.name
    async with result.log.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    assert "Reading data... done." in text, "test.log doesn't report reading the data"


@env.task
async def test_humid_paired_end() -> None:
    # upstream case: humid "sarscov2 - paired-end"
    trimmed = await fastp(
        await fixture(FASTQ + "test_1.fastq.gz"),
        reads_2=await fixture(FASTQ + "test_2.fastq.gz"),
        prefix="test",
        args=FASTP_ARGS,
    )
    result = await humid(trimmed.reads[0], reads_2=trimmed.reads[1], prefix="test", args=HUMID_ARGS)
    await check_gunzipped(
        result.dedup,
        {
            "test_R1_dedup.fastp.fastq.gz": "9c736ed873cf976039886229b742e272",
            "test_R2_dedup.fastp.fastq.gz": "730c5450425a26230a5377972b84ba4c",
        },
        "dedup",
    )
    await check_gunzipped(
        result.annotated,
        {
            "test_R1_annotated.fastp.fastq.gz": "24bdbebb3ae1b3e96905df28c61e8571",
            "test_R2_annotated.fastp.fastq.gz": "2bcacbd846fb74e2f19985967ce30c0a",
        },
        "annotated",
    )
    await check_stats(
        result.stats,
        {
            "clusters.dat": "b88d9d3330fa58f35800ce1df9755509",
            "counts.dat": "b88d9d3330fa58f35800ce1df9755509",
            "neigh.dat": "c4780d93e68a773dc1cf29b50319c476",
            "stats.dat": "1bcd2377bc9261e8466123af79b74c33",
        },
    )
    await check_log(result)


@env.task
async def test_humid_single_end() -> None:
    # upstream case: humid "sarscov2 - single-end"
    trimmed = await fastp(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", args=FASTP_ARGS)
    result = await humid(trimmed.reads[0], prefix="test", args=HUMID_ARGS)
    await check_gunzipped(result.dedup, {"test_dedup.fastp.fastq.gz": "44ff29f56dd5d48e376663bc22fbe3b3"}, "dedup")
    await check_gunzipped(
        result.annotated, {"test_annotated.fastp.fastq.gz": "567bf1de24726825a0f9a9a0aede60d2"}, "annotated"
    )
    await check_stats(
        result.stats,
        {
            "clusters.dat": "6b559f508be032a6f43a6135776655c5",
            "counts.dat": "6b559f508be032a6f43a6135776655c5",
            "neigh.dat": "a2e34a81a3cbda2cd61005fa66c31961",
            "stats.dat": "17d662b83217a6537247903940271bef",
        },
    )
    await check_log(result)


tests = [test_humid_paired_end, test_humid_single_end]
