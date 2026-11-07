"""Tests for flyte_bio.modules.sam2lca.

The case is upstream's sam2lca/analyze test (modules commit 9a1cc9f, sam2lca
1.1.4): sam2lca's own test alignment against its ``test`` taxonomy and
accession-to-taxon databases (``-a test -t test``), which sam2lca builds into
an empty database directory, downloading them. Upstream's snapshot marks the
CSV and JSON unstable and checks only their names (and that no BAM is
written); here the test also checks they hold the same, non-empty set of taxa.
"""

import csv
import io
import json

from flyte.io import File

from flyte_bio.modules.sam2lca import sam2lca_analyze
from tests.framework import env, fixture

MICROTEST = "https://raw.githubusercontent.com/maxibor/sam2lca/1.1.2/tests/data/microtest.sorted.bam"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_sam2lca_analyze() -> None:
    # upstream case: "test-sam2lca-analyze"
    result = await sam2lca_analyze(
        await fixture(MICROTEST), await fixture(MICROTEST + ".bai"), prefix="test", args="-a test -t test"
    )
    assert result.csv.path.endswith("/test.sam2lca.csv"), f"csv is {result.csv.path}"
    assert result.json.path.endswith("/test.sam2lca.json"), f"json is {result.json.path}"
    assert result.bam == [], f"sam2lca wrote BAMs without -b: {[f.path for f in result.bam]}"

    rows = list(csv.DictReader(io.StringIO(await text(result.csv))))
    assert rows, "test.sam2lca.csv has no taxa"
    assert "TAXID" in rows[0] and "count_descendant" in rows[0], f"unexpected csv columns {list(rows[0])}"
    taxa = json.loads(await text(result.json))
    assert sorted(taxa) == sorted(r["TAXID"] for r in rows), "the csv and json list different taxa"


tests = [test_sam2lca_analyze]
