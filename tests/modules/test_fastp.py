"""Tests for flyte_bio.modules.fastp.

The cases and expected md5s are upstream's fastp tests (modules commit 9a1cc9f,
fastp 1.4.0). nf-test computes the
md5s over the decompressed reads. The module's 4 CPUs match the 4 that
upstream's thread-checking cases set.
"""

from flyte.io import File

from flyte_bio.modules.fastp import FastpResult, fastp
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
LITTLE_ADAPTER = "The input has little adapter percentage (~0.000000%), probably it's trimmed before."


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check_reports(result: FastpResult, html: str, log: str, json: str | None = None) -> None:
    assert html in await text(result.html), f"fastp html lacks {html!r}"
    assert log in await text(result.log), f"fastp log lacks {log!r}"
    if json is not None:
        assert json in await text(result.json), f"fastp json lacks {json!r}"


async def check_reads(files: list[File], md5s: list[str], label: str) -> None:
    assert len(files) == len(md5s), f"{label}: {len(files)} files, expected {len(md5s)}"
    for file, md5 in zip(files, md5s, strict=True):
        await assert_gunzipped_md5(file, md5, label=f"{label} {file.name}")


async def paired() -> tuple[File, File]:
    return await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")


@env.task
async def test_fastp_single_end() -> None:
    result = await fastp(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    await check_reads(result.reads, ["67b2bbae47f073e05a97a9c2edce23c7"], "single-end")
    assert result.reads_fail == [] and result.reads_merged is None
    await check_reports(result, "single end (151 cycles)", "reads passed filter: 99", "--thread 3 ")


@env.task
async def test_fastp_paired_end() -> None:
    result = await fastp(*await paired(), prefix="test")
    await check_reads(
        result.reads, ["67b2bbae47f073e05a97a9c2edce23c7", "25cbdca08e2083dbd4f0502de6b62f39"], "paired-end"
    )
    await check_reports(result, LITTLE_ADAPTER, "Q30 bases: 12281(88.3716%)", "--thread 2 ")


@env.task
async def test_fastp_interleaved() -> None:
    reads = await fixture(FASTQ + "test_interleaved.fastq.gz")
    result = await fastp(reads, prefix="test", args="--interleaved_in -e 30")
    await check_reads(result.reads, ["217d62dc13a23e92513a1bd8e1bcea39"], "interleaved")
    await check_reports(result, "paired end (151 cycles + 151 cycles)", "reads passed filter: 162", "--thread 2 ")


@env.task
async def test_fastp_single_end_trim_fail() -> None:
    result = await fastp(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", save_trimmed_fail=True)
    await check_reads(result.reads, ["67b2bbae47f073e05a97a9c2edce23c7"], "single-end")
    await check_reads(result.reads_fail, ["3e4aaadb66a5b8fc9b881bf39c227abd"], "single-end failed")
    await check_reports(result, "single end (151 cycles)", "reads passed filter: 99")


@env.task
async def test_fastp_paired_end_trim_fail() -> None:
    result = await fastp(*await paired(), prefix="test", save_trimmed_fail=True, args="-e 30")
    await check_reads(
        result.reads, ["6ff32a64c5188b9a9192be1398c262c7", "db0cb7c9977e94ac2b4b446ebd017a8a"], "paired-end"
    )
    await check_reads(
        result.reads_fail,
        ["409b687c734cedd7a1fec14d316e1366", "4f273cf3159c13f79e8ffae12f5661f6", "f97b9edefb5649aab661fbc9e71fc995"],
        "paired-end failed",
    )
    await check_reports(result, LITTLE_ADAPTER, "reads passed filter: 162")


@env.task
async def test_fastp_paired_end_merged() -> None:
    result = await fastp(*await paired(), prefix="test", save_merged=True)
    await check_reads(
        result.reads, ["54b726a55e992a869fd3fa778afe1672", "29d3b33b869f7b63417b8ff07bb128ba"], "paired-end"
    )
    assert result.reads_merged is not None
    await check_reads([result.reads_merged], ["c873bb1ab3fa859dcc47306465e749d5"], "merged")
    await check_reports(result, LITTLE_ADAPTER, "total reads: 75")


@env.task
async def test_fastp_paired_end_merged_adapterlist() -> None:
    adapters = await fixture("delete_me/fastp/adapters.fasta")
    result = await fastp(*await paired(), adapter_fasta=adapters, prefix="test", save_merged=True)
    await check_reads(
        result.reads, ["54b726a55e992a869fd3fa778afe1672", "29d3b33b869f7b63417b8ff07bb128ba"], "paired-end"
    )
    assert result.reads_merged is not None
    await check_reads([result.reads_merged], ["c873bb1ab3fa859dcc47306465e749d5"], "merged")
    await check_reports(result, "<div id='After_filtering__merged__quality'>", "total bases: 13683")


@env.task
async def test_fastp_single_end_qc_only() -> None:
    result = await fastp(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", discard_trimmed_pass=True)
    assert result.reads == [] and result.reads_fail == [] and result.reads_merged is None
    assert await result.results.get_file("test.fastp.fastq.gz") is None, "qc-only run wrote reads"
    await check_reports(result, "single end (151 cycles)", "reads passed filter: 99")


@env.task
async def test_fastp_paired_end_qc_only() -> None:
    result = await fastp(*await paired(), prefix="test", discard_trimmed_pass=True)
    assert result.reads == []
    assert await result.results.get_file("test_R1.fastp.fastq.gz") is None, "qc-only run wrote reads"
    await check_reports(result, LITTLE_ADAPTER, "Q30 bases: 12281(88.3716%)")


tests = [
    test_fastp_single_end,
    test_fastp_paired_end,
    test_fastp_interleaved,
    test_fastp_single_end_trim_fail,
    test_fastp_paired_end_trim_fail,
    test_fastp_paired_end_merged,
    test_fastp_paired_end_merged_adapterlist,
    test_fastp_single_end_qc_only,
    test_fastp_paired_end_qc_only,
]
