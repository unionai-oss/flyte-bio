"""Tests for flyte_bio.modules.krona.

The case is upstream's krona/ktimporttaxonomy test (modules commit 9a1cc9f,
KronaTools 2.8.1). The chart embeds KronaTools' resource URLs and isn't
snapshotted upstream, which checks only its name; this test also checks it's a
Krona chart.

The ktUpdateTaxonomy case is upstream's krona/ktupdatetaxonomy test (same
commit). NCBI's taxonomy changes daily, so upstream checks only the
database's name; this test also checks it's a KronaTools taxonomy, rooted at
taxon 1.
"""

from flyte_bio.modules.krona import krona_ktimporttaxonomy, krona_ktupdatetaxonomy
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


@env.task
async def test_krona_ktupdatetaxonomy() -> None:
    taxonomy = await krona_ktupdatetaxonomy()
    assert taxonomy.name == "taxonomy.tab", taxonomy.name
    async with taxonomy.open("rb") as fh:
        head = bytes(await fh.read(1 << 16)).decode(errors="replace")
    root = head.splitlines()[0].split("\t")
    # taxonomy.tab rows: taxon ID, depth, parent, rank, name.
    assert root[0] == "1" and root[-1] == "root", f"taxonomy.tab starts with {root}"


tests = [test_krona_ktimporttaxonomy_kraken_report, test_krona_ktupdatetaxonomy]
