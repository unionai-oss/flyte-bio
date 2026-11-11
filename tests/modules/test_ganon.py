"""Tests for flyte_bio.modules.ganon.

The cases and expected md5s are upstream's ganon/buildcustom tests (modules
commit 9a1cc9f, ganon 2.1.0), with each test's ``--input-target``. In the
FASTA case ganon downloads the NCBI taxonomy, so, as upstream, its
``test.tax`` is only checked for the genome's accession.

The classify, report and table cases are upstream's ganon/classify,
ganon/report and ganon/table tests (same commit), with their configs'
options, over the FASTA case's database (built with ``--input-target
sequence``, as their setups).
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.ganon import (
    GanonBuildResult,
    GanonClassifyResult,
    ganon_buildcustom,
    ganon_classify,
    ganon_report,
    ganon_table,
)
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"
METAGENOME = "genomics/sarscov2/metagenome/"


async def ganon_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_ganon_buildcustom_fasta() -> None:
    # upstream case: ganon/buildcustom "sarscov2 - genome - fasta"
    result = await ganon_buildcustom(await fixture(GENOME), prefix="test", args="--input-target sequence")
    assert result.filter.name == "test.hibf", result.filter.name
    await assert_md5(result.filter, "d10fe6fc6d198696bc15ca85a1459614", label="test.hibf")
    assert result.tax is not None, "ganon wrote no test.tax"
    assert "MT192765.1" in await ganon_text(result.tax), "test.tax lacks MT192765.1"


@env.task
async def test_ganon_buildcustom_tsv() -> None:
    # upstream case: ganon/buildcustom "sarscov2 - genome - tsv"
    tsv = Path(tempfile.mkdtemp()) / "ganon_fasta_input.tsv"
    tsv.write_text("genome.fasta\tSevere_acute_respiratory_syndrome_coronavirus_2\t2697049\n")
    result = await ganon_buildcustom(
        await fixture(GENOME),
        input_tsv=await File.from_local(str(tsv)),
        taxonomy_files=[await fixture(METAGENOME + "prot_nodes.dmp"), await fixture(METAGENOME + "prot_names.dmp")],
        prefix="test",
        args="--input-target file",
    )
    await assert_md5(result.filter, "9edfe4c3873d621a88ebcad438dca42c", label="test.hibf")
    assert result.tax is not None, "ganon wrote no test.tax"
    await assert_md5(result.tax, "5fd3dd6f618d01b33d2ad1c3884aba1d", label="test.tax")


FASTQ = "genomics/sarscov2/illumina/fastq/"
CLASSIFY_ARGS = "--output-one --output-all --output-unclassified"


async def ganon_test_db() -> GanonBuildResult:
    return await ganon_buildcustom(await fixture(GENOME), prefix="test", args="--input-target sequence")


async def ganon_test_classify(db: GanonBuildResult, paired: bool) -> GanonClassifyResult:
    reads_2 = await fixture(FASTQ + "test_2.fastq.gz") if paired else None
    return await ganon_classify(
        await fixture(FASTQ + "test_1.fastq.gz"), db.db, reads_2=reads_2, prefix="test", args=CLASSIFY_ARGS
    )


async def check_ganon_classify(result: GanonClassifyResult, tre_md5: str, rep_md5: str) -> None:
    await assert_md5(result.tre, tre_md5, label="test.tre")
    await assert_md5(result.rep, rep_md5, label="test.rep")
    assert "Total elapsed time:" in await ganon_text(result.log), "test.log lacks 'Total elapsed time:'"


@env.task
async def test_ganon_classify_single_end() -> None:
    # upstream case: ganon/classify "sarscov2 single-end [fastq]"
    result = await ganon_test_classify(await ganon_test_db(), paired=False)
    await check_ganon_classify(result, "828f3ea0e15ae28d95f71ae4874edc38", "9b5f8051aa459a80a678e36259a15746")


@env.task
async def test_ganon_classify_paired_end() -> None:
    # upstream case: ganon/classify "sarscov2 paired-end [fastq]"
    result = await ganon_test_classify(await ganon_test_db(), paired=True)
    await check_ganon_classify(result, "2eebb02e0e37268314d41cc274ec3aa0", "757e5cbf7a978913277401285d743336")


REPORT_ARGS = "--report-type reads --output-format tsv"


@env.task
async def test_ganon_report_single_end() -> None:
    # upstream case: ganon/report "sarscov2 single-end [fastq]"
    db = await ganon_test_db()
    classified = await ganon_test_classify(db, paired=False)
    tre = await ganon_report(classified.rep, db.db, prefix="test", args=REPORT_ARGS)
    assert tre.name == "test.tre", tre.name
    await assert_md5(tre, "828f3ea0e15ae28d95f71ae4874edc38", label="test.tre")


@env.task
async def test_ganon_table_single_end() -> None:
    # upstream case: ganon/table "sarscov2 single-end [fastq]"
    db = await ganon_test_db()
    classified = await ganon_test_classify(db, paired=False)
    tre = await ganon_report(classified.rep, db.db, prefix="test", args=REPORT_ARGS)
    table = await ganon_table([tre], prefix="db1", args="--output-format tsv --output-value percentages")
    assert table.name == "db1.txt", table.name
    await assert_md5(table, "1da6f44d95981cca78a50b32d19f2337", label="db1.txt")


tests = [
    test_ganon_buildcustom_fasta,
    test_ganon_buildcustom_tsv,
    test_ganon_classify_single_end,
    test_ganon_classify_paired_end,
    test_ganon_report_single_end,
    test_ganon_table_single_end,
]
