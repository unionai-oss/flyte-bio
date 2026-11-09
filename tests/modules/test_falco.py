"""Tests for flyte_bio.modules.falco.

The cases are upstream's falco tests (modules commit 9a1cc9f, falco 1.2.5).
Upstream snapshots only the reports' names, so these check the names, and
that each data table is a falco report (``##Falco`` header).
"""

import os

from flyte.io import File

from flyte_bio.modules.falco import FalcoResult, falco
from tests.framework import env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


def falco_names(files: list[File]) -> list[str]:
    return [os.path.basename(f.path) for f in files]


async def check_falco(result: FalcoResult, txt: list[str], html: list[str]) -> None:
    assert falco_names(result.txt) == txt, f"txt reports {falco_names(result.txt)}"
    assert falco_names(result.html) == html, f"html reports {falco_names(result.html)}"
    for table in result.txt:
        if table.path.endswith("_fastqc_data.txt"):
            async with table.open("rb") as fh:
                head = bytes(await fh.read(7))
            assert head == b"##Falco", f"{table.path} starts with {head!r}"


@env.task
async def test_falco_single_end() -> None:
    # upstream case: "sarscov2 - fastq - single end"
    result = await falco([await fixture(FASTQ + "test_1.fastq.gz")], prefix="test")
    await check_falco(result, ["test_fastqc_data.txt", "test_summary.txt"], ["test_report.html"])


@env.task
async def test_falco_paired_end() -> None:
    # upstream case: "sarscov2 - fastq - paired end"
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    result = await falco(reads, prefix="test")
    await check_falco(
        result,
        [
            "test_1.fastq.gz_fastqc_data.txt",
            "test_1.fastq.gz_summary.txt",
            "test_2.fastq.gz_fastqc_data.txt",
            "test_2.fastq.gz_summary.txt",
        ],
        ["test_1.fastq.gz_fastqc_report.html", "test_2.fastq.gz_fastqc_report.html"],
    )


@env.task
async def test_falco_interleaved() -> None:
    # upstream case: "sarscov2 - fastq - interleaved"
    result = await falco([await fixture(FASTQ + "test_interleaved.fastq.gz")], prefix="test")
    await check_falco(result, ["test_fastqc_data.txt", "test_summary.txt"], ["test_report.html"])


tests = [test_falco_single_end, test_falco_paired_end, test_falco_interleaved]
