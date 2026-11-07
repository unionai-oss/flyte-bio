"""Tests for flyte_bio.modules.krona.

The case is upstream's krona/ktimporttaxonomy test (modules commit 9a1cc9f,
KronaTools 2.8.1). The chart embeds KronaTools' resource URLs and isn't
snapshotted upstream, which checks only its name; this test also checks it's a
Krona chart.
"""

from flyte_bio.modules.krona import krona_ktimporttaxonomy
from tests.framework import env, fixture

METAGENOME = "genomics/sarscov2/metagenome/"


@env.task
async def test_krona_ktimporttaxonomy_kraken_report() -> None:
    report = await fixture(METAGENOME + "test_1.kraken2.report.txt")
    taxonomy = await fixture(METAGENOME + "krona_taxonomy.tab")
    html = await krona_ktimporttaxonomy(report, taxonomy, prefix="test")
    assert html.name == "test.html", html.name
    async with html.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    assert "<krona" in text, "test.html isn't a Krona chart"


tests = [test_krona_ktimporttaxonomy_kraken_report]
