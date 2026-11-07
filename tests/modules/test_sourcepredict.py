"""Tests for flyte_bio.modules.sourcepredict.

The case is upstream's sourcepredict test (modules commit 9a1cc9f,
SourcePredict 0.5.1, ``-me MDS``), whose setup decompresses the taxonomy
database with xz/decompress, as here. Upstream doesn't pin the predicted
proportions (floating-point differences) and its snapshot of the table's
structure is empty; here the test checks the table's shape instead: a row per
source label plus ``unknown``, a column for the sink sample, and proportions
that sum to 1.
"""

import csv
import io

from flyte_bio.modules.sourcepredict import sourcepredict
from flyte_bio.modules.xz import xz_decompress
from tests.framework import env, fixture

MISC = "genomics/prokaryotes/metagenome/taxonomy/misc/"


@env.task
async def test_sourcepredict_metagenome_csv() -> None:
    # upstream case: "metagenome - csv"
    taxa_sqlite = await xz_decompress(await fixture(MISC + "taxa_sqlite.xz"))
    report = await sourcepredict(
        await fixture(MISC + "sink_taxid.csv"),
        sources=await fixture(MISC + "sources_taxid.csv"),
        labels=await fixture(MISC + "sources_labels.csv"),
        taxa_sqlite=taxa_sqlite,
        taxa_sqlite_traverse_pkl=await fixture(MISC + "taxa_sqlite_traverse.pkl"),
        prefix="test",
        args="-me MDS",
    )
    assert report.path.endswith("/test.sourcepredict.csv"), f"report is {report.path}"
    async with report.open("rb") as fh:
        rows = list(csv.reader(io.StringIO(bytes(await fh.read()).decode())))
    assert rows[0][1:] == ["ERR1915662"], f"unexpected sink columns {rows[0]}"
    sources = sorted(r[0] for r in rows[1:])
    assert sources == ["Canis_familiaris", "Homo_sapiens", "unknown"], f"unexpected sources {sources}"
    total = sum(float(r[1]) for r in rows[1:])
    assert abs(total - 1) < 1e-6, f"the proportions sum to {total}, not 1"


tests = [test_sourcepredict_metagenome_csv]
