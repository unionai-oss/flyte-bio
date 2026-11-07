"""Tests for flyte_bio.modules.wipertools.

The cases and expected values are upstream's wipertools fastqscatter,
fastqwiper, fastqgather and reportgather tests (modules commit 9a1cc9f,
wipertools 1.1.4), with their test configs' ``args``. nf-test hashes
``.fastq.gz`` outputs decompressed.

Upstream's reportgather case mixes the reports of two fastqwiper runs into a
channel without collecting it, so reportgather runs once per report; these
tests do the same, one gather per report. Upstream checks the line count,
header and well-formed-lines line of the first gathered report, which is
test_truncated_clean.fastq's.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.wipertools import (
    wipertools_fastqgather,
    wipertools_fastqscatter,
    wipertools_fastqwiper,
    wipertools_reportgather,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

SARSCOV2_FASTQ = "genomics/sarscov2/illumina/fastq/test_1.fastq.gz"
QUALITY_MISMATCH = "genomics/homo_sapiens/illumina/fastq/test_quality_mismatch.fastq"
TRUNCATED_CLEAN = "genomics/homo_sapiens/illumina/fastq/test_truncated_clean.fastq"
QUALITY_MISMATCH_REPORT_MD5 = "d6fe1f0066f538993bb8f6b90150c8d5"


async def lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        data = bytes(await fh.read())
    if file.name.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode().splitlines()


@env.task
async def test_wipertools_fastqscatter_fastq() -> None:
    # upstream case: wipertools/fastqscatter "test1 - fastq"
    result = await wipertools_fastqscatter(
        await fixture(SARSCOV2_FASTQ), 2, prefix="test", args="-o chunks_custom -O unix"
    )
    names = [f.name for f in result.chunks]
    assert names == ["test_000.fastq.gz", "test_001.fastq.gz"], names
    await assert_gunzipped_md5(result.chunks[0], "ecc4b1841cd94704bba742ea4dcd48b0", label=names[0])
    await assert_gunzipped_md5(result.chunks[1], "b3de467f2b6ab0d14e1f6ce14932a411", label=names[1])


@env.task
async def test_wipertools_fastqwiper_quality_mismatch() -> None:
    # upstream case: wipertools/fastqwiper "Quality mismatch - fastq"
    result = await wipertools_fastqwiper(await fixture(QUALITY_MISMATCH), prefix="test_wiped", args="-a ACGT")
    assert result.fastq.name == "test_wiped.fastq.gz", result.fastq.name
    reads = await lines(result.fastq)
    assert len(reads) == 8, reads
    assert reads[7] == "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", reads[7]
    report = await lines(result.report)
    assert len(report) == 11, report
    assert report[0] == "FASTQWIPER REPORT:", report[0]
    assert report[9] == "Reads discarded because len(SEQ) neq len(QUAL): 1", report[9]
    await assert_gunzipped_md5(result.fastq, "81b1d99cf5c1615123a8c31760858a10", label="wiped reads")
    await assert_md5(result.report, QUALITY_MISMATCH_REPORT_MD5, label="report")


@env.task
async def test_wipertools_fastqgather_two_chunks() -> None:
    # upstream case: wipertools/fastqgather "merge two gzipped fastq - fastq.gz"
    scattered = await wipertools_fastqscatter(await fixture(SARSCOV2_FASTQ), 2, prefix="fastq_chunks")
    gathered = await wipertools_fastqgather(
        scattered.chunks, prefix="fastq_chunks_gather", args="--prefix fastq_chunks"
    )
    assert gathered.name == "fastq_chunks_gather.fastq.gz", gathered.name
    await assert_gunzipped_md5(gathered, "4161df271f9bfcd25d5845a1e220dbec", label="gathered reads")


@env.task
async def test_wipertools_reportgather_two_reports() -> None:
    # upstream case: wipertools/reportgather "merge two reports - .report" (one gather per report)
    mismatch = await wipertools_fastqwiper(await fixture(QUALITY_MISMATCH), prefix="test1_wiped")
    truncated = await wipertools_fastqwiper(await fixture(TRUNCATED_CLEAN), prefix="test2_wiped")
    first = await wipertools_reportgather([truncated.report], prefix="reports_gather")
    second = await wipertools_reportgather([mismatch.report], prefix="reports_gather")
    assert first.name == "reports_gather.report", first.name
    report = await lines(first)
    assert len(report) == 11, report
    assert report[0] == "FASTQWIPER REPORT:", report[0]
    assert report[3] == "Well-formed lines: 8 (80.0%)", report[3]
    await assert_md5(first, "3d2cfb8d7b72afe9b88663cf9dbf8129", label="truncated-clean report")
    await assert_md5(second, QUALITY_MISMATCH_REPORT_MD5, label="quality-mismatch report")


tests = [
    test_wipertools_fastqscatter_fastq,
    test_wipertools_fastqwiper_quality_mismatch,
    test_wipertools_fastqgather_two_chunks,
    test_wipertools_reportgather_two_reports,
]
