"""Tests for flyte_bio.modules.barrnap.

The case and expected md5 are upstream's barrnap test (modules commit
9a1cc9f, Barrnap 0.9).
"""

from flyte_bio.modules.barrnap import barrnap
from tests.framework import assert_md5, env, fixture


@env.task
async def test_barrnap_contigs() -> None:
    gff = await barrnap(await fixture("genomics/sarscov2/illumina/fasta/contigs.fasta"), kingdom="bac", prefix="test")
    await assert_md5(gff, "df19e1b84ba6f691d20c72b397c88abf", label="gff")


tests = [test_barrnap_contigs]
