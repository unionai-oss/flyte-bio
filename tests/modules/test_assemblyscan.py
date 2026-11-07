"""Tests for flyte_bio.modules.assemblyscan.

The cases and expected values are upstream's assemblyscan tests (modules
commit 9a1cc9f, assembly-scan 1.0.0): the report of the SARS-CoV-2 genome as
TSV, and as JSON with ``--json`` (upstream's ``json.config``).
"""

from flyte_bio.modules.assemblyscan import assemblyscan
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"


@env.task
async def test_assemblyscan_tsv() -> None:
    # upstream case: assemblyscan "genome.fasta"
    report = await assemblyscan(await fixture(GENOME), prefix="test")
    assert report.path.endswith("/test.tsv"), report.path
    await assert_md5(report, "d3bda95a4a6f65cd94adeccca998b0b6", label="test.tsv")


@env.task
async def test_assemblyscan_json() -> None:
    # upstream case: assemblyscan "genome.fasta - json"
    report = await assemblyscan(await fixture(GENOME), prefix="test", args="--json")
    assert report.path.endswith("/test.json"), report.path
    await assert_md5(report, "bcef4e2bf296339c1ca7c513b483de96", label="test.json")


tests = [test_assemblyscan_tsv, test_assemblyscan_json]
