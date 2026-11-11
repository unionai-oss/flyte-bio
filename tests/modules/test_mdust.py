"""Tests for flyte_bio.modules.mdust.

The case and md5 are upstream's mdust test (modules commit 9a1cc9f, mdust
2006.10.17). The SARS-CoV-2 genome has no low-complexity stretch, so the
output equals the input.
"""

from flyte_bio.modules.mdust import mdust
from tests.framework import assert_md5, env, fixture


@env.task
async def test_mdust_sarscov2() -> None:
    # upstream case: mdust "sarscov2 - fasta"
    masked = await mdust(await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test_mdust")
    await assert_md5(masked, "483f4a5dfe60171c86ee9b7e6dff908b", label="test_mdust.fasta")


tests = [test_mdust_sarscov2]
