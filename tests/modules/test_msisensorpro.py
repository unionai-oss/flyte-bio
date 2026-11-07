"""Tests for flyte_bio.modules.msisensorpro.

The case and its expected md5 are upstream's msisensorpro/scan test (modules
commit 9a1cc9f, MSIsensor-pro 1.3.0).
"""

from flyte_bio.modules.msisensorpro import msisensorpro_scan
from tests.framework import assert_md5, env, fixture


@env.task
async def test_msisensorpro_scan_fasta() -> None:
    # upstream case: msisensorpro/scan "test-msisensorpro-scan"
    fasta = await fixture("genomics/homo_sapiens/genome/chr21/sequence/genome.fasta")
    sites = await msisensorpro_scan(fasta, prefix="test")
    assert sites.name == "test.msisensor_scan.list", sites.name
    await assert_md5(sites, "614754c7f1f44d5988dd80f8f21f69d1", label="scan list")


tests = [test_msisensorpro_scan_fasta]
