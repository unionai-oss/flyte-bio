"""Tests for flyte_bio.modules.fastavalidator.

The cases and expected md5s are upstream's fastavalidator tests at modules
commit 41dfa3f, the module's last working version (py_fasta_validator 0.6). At
9a1cc9f upstream's only test checks that the deprecated module fails.
"""

from flyte_bio.modules.fastavalidator import fastavalidator
from tests.framework import assert_md5, env, fixture


async def read_text(file) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_fastavalidator_sarscov2_fasta_valid() -> None:
    result = await fastavalidator(await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test")
    assert result.error_log is None, "a valid FASTA wrote an error log"
    assert result.success_log is not None, "no success log"
    assert "Validation successful..." in await read_text(result.success_log)
    await assert_md5(result.success_log, "b0b859eda1db5cd43915846e00ebc22c", label="test.success.log")


@env.task
async def test_fastavalidator_sarscov2_gff3_invalid() -> None:
    result = await fastavalidator(await fixture("genomics/sarscov2/genome/genome.gff3"), prefix="test")
    assert result.success_log is None, "an invalid FASTA wrote a success log"
    assert result.error_log is not None, "no error log"
    assert "genome.gff3 does not start with a >" in await read_text(result.error_log)
    await assert_md5(result.error_log, "531d520c0e7767176f743f197f1f87b3", label="test.error.log")


tests = [test_fastavalidator_sarscov2_fasta_valid, test_fastavalidator_sarscov2_gff3_invalid]
