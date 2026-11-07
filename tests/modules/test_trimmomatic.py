"""Tests for flyte_bio.modules.trimmomatic.

The cases are upstream's trimmomatic tests (modules commit 9a1cc9f,
Trimmomatic 0.39), with the steps their configs set and their snapshot md5s
(over the decompressed reads).
"""

from flyte_bio.modules.trimmomatic import TrimmomaticResult, trimmomatic
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
STEPS = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"


async def check(result: TrimmomaticResult, reads: dict[str, str], summary: str, label: str) -> None:
    assert [f.name for f in result.reads] == list(reads), f"{label}: reads {[f.name for f in result.reads]}"
    for f, md5 in zip(result.reads, reads.values(), strict=True):
        await assert_gunzipped_md5(f, md5, label=f"{label} {f.name}")
    assert result.summary.name == "test.summary", result.summary.name
    await assert_md5(result.summary, summary, label=f"{label} summary")
    assert result.trim_log.name == "test_trim.log", result.trim_log.name


@env.task
async def test_trimmomatic_single_end() -> None:
    # upstream case: trimmomatic "Single-Read" (nextflow_SE.config)
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await trimmomatic(reads, prefix="test", args=f"ILLUMINACLIP:TruSeq3-SE:2:30:10 {STEPS}")
    assert result.unpaired_reads == []
    await check(
        result,
        {"test.SE.paired.trim.fastq.gz": "e68abbd3b88f7ec12940a4f5c2b8bfb9"},
        "24c973237557a1439c775ca19a5deaa5",
        "single-end",
    )


@env.task
async def test_trimmomatic_paired_end() -> None:
    # upstream case: trimmomatic "Paired-Reads" (nextflow_PE.config)
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    result = await trimmomatic(r1, r2, prefix="test", args=f"ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 {STEPS}")
    assert [f.name for f in result.unpaired_reads] == ["test.unpaired.trim_1.fastq.gz", "test.unpaired.trim_2.fastq.gz"]
    await check(
        result,
        {
            "test.paired.trim_1.fastq.gz": "a866e64f451745f176316d0df7d52b30",
            "test.paired.trim_2.fastq.gz": "725d4ab909b39291ab56b090cab32075",
        },
        "9698e5e5c060bbe64588998fe35f8d71",
        "paired-end",
    )


@env.task
async def test_trimmomatic_no_steps_fails() -> None:
    # upstream case: trimmomatic "No Adaptors": without trimming steps the process fails
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    try:
        await trimmomatic(r1, r2, prefix="test")
    except Exception:
        return
    raise AssertionError("No Adaptors: Trimmomatic ran without any trimming steps")


tests = [test_trimmomatic_single_end, test_trimmomatic_paired_end, test_trimmomatic_no_steps_fails]
