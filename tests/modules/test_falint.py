"""Tests for flyte_bio.modules.falint.

The cases and expected md5 are upstream's falint tests (modules commit
9a1cc9f, fa-lint 1.2.0): a valid FASTA, and a GFF3 that fa-lint rejects.
"""

from flyte.io import File

from flyte_bio.modules.falint import falint
from tests.framework import assert_md5, env, fixture


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_falint_fasta_valid() -> None:
    result = await falint(await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test")
    assert result.error_log is None, "fa-lint wrote an error log for a valid FASTA"
    assert result.success_log is not None, "fa-lint wrote no success log"
    await assert_md5(result.success_log, "7b9802635c4f44d1711e7132ccdb648e", label="success log")


@env.task
async def test_falint_gff3_invalid() -> None:
    result = await falint(await fixture("genomics/sarscov2/genome/genome.gff3"), prefix="test")
    assert result.success_log is None, "fa-lint wrote a success log for a GFF3"
    assert result.error_log is not None, "fa-lint wrote no error log"
    message = "Fasta file must start with a header line"
    assert message in await text(result.error_log), f"error log lacks {message!r}"


tests = [test_falint_fasta_valid, test_falint_gff3_invalid]
