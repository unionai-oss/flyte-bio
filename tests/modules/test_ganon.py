"""Tests for flyte_bio.modules.ganon.

The cases and expected md5s are upstream's ganon/buildcustom tests (modules
commit 9a1cc9f, ganon 2.1.0), with each test's ``--input-target``. In the
FASTA case ganon downloads the NCBI taxonomy, so, as upstream, its
``test.tax`` is only checked for the genome's accession.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.ganon import ganon_buildcustom
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


tests = [test_ganon_buildcustom_fasta, test_ganon_buildcustom_tsv]
