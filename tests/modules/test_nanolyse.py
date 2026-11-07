"""Tests for flyte_bio.modules.nanolyse.

The case and its expected md5 are upstream's nanolyse test (modules commit
9a1cc9f, NanoLyse 1.2.0), with its test config's ``<id>.clean`` prefix.
nf-test hashes the decompressed reads; the log, which records timestamps, is
snapshotted by name only.
"""

from flyte_bio.modules.nanolyse import nanolyse
from tests.framework import assert_gunzipped_md5, env, fixture


@env.task
async def test_nanolyse_lambda() -> None:
    # upstream case: nanolyse "test-nanolyse"
    fastq = await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz")
    lambda_fasta = await fixture("https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz")
    result = await nanolyse(fastq, lambda_fasta, prefix="test.clean")
    assert result.reads.name == "test.clean.fastq.gz", result.reads.name
    assert result.log.name == "test.clean.nanolyse.log", result.log.name
    await assert_gunzipped_md5(result.reads, "7567d853ada6ac142332619d0b541d76", label="clean reads")
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    assert "NanoLyse: removed 0 reads." in log, log


tests = [test_nanolyse_lambda]
