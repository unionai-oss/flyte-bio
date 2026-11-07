"""Tests for flyte_bio.modules.qcatch.

Upstream's test data: a simpleaf quantification of 1k PBMCs (10x v3) from
QCatch's authors. Upstream snapshots the md5 of the filtered h5ad and of the
metrics summary; the doublet run only checks that QCatch succeeds.

QCatch's cell calling isn't deterministic: two identical runs kept 1,226 and
1,229 cells (upstream's metrics md5 is 22deb2f3aea4e4abf476970576362428), so
the summary is checked by content: its inputs exactly, the retained cells and
their per-cell medians within narrow bounds.
"""

import csv
import io

from flyte.io import File

from flyte_bio.modules.qcatch import qcatch
from tests.framework import assert_nonempty, env, fixture_dir

SUMMARY_HEADER = [
    "Number of retained cells",
    "Number of all processed cells",
    "Mean reads per retained cell",
    "Median UMI per retained cell",
    "Median genes per retained cell",
    "Total genes detected for retained cells",
    "Mapping rate",
    "Sequencing saturation",
]


async def check_metrics_summary(summary: File) -> None:
    async with summary.open("rb") as fh:
        rows = list(csv.reader(io.StringIO(bytes(await fh.read()).decode())))
    assert rows[0] == SUMMARY_HEADER, rows[0]
    values = dict(zip(rows[0], rows[1], strict=True))

    def number(name: str) -> int:
        return int(values[name].replace(",", ""))

    # What QCatch reads from the quantification is fixed.
    assert number("Number of all processed cells") == 76077, values
    assert values["Mapping rate"] == "N/A", values
    assert values["Sequencing saturation"] == "70.23%", values
    # Cell calling varies slightly between runs (1,226 and 1,229 seen).
    assert 1200 <= number("Number of retained cells") <= 1260, values
    for name, expected in (
        ("Mean reads per retained cell", 46758),
        ("Median UMI per retained cell", 10486),
        ("Median genes per retained cell", 3276),
        ("Total genes detected for retained cells", 26374),
    ):
        assert abs(number(name) - expected) <= expected * 0.02, (name, values)


TEST_DATA = "https://umd.box.com/shared/static/zd4sai70uw9fs24e1qx6r41ec50pf45g.zip"


@env.task
async def test_qcatch() -> None:
    quant = await fixture_dir(TEST_DATA, "test_data/1k_pbmc_v3")
    result = await qcatch(quant, prefix="1k_pbmc_v3", chemistry="10X_3p_v3")
    assert result.report.path.endswith("/1k_pbmc_v3_qcatch_report.html"), result.report.path
    await assert_nonempty(result.report, label="report")
    # The h5ad's bytes vary between runs (HDF5), so unlike upstream only its existence is checked.
    await assert_nonempty(result.filtered_h5ad, label="filtered h5ad")
    await check_metrics_summary(result.metrics_summary)


@env.task
async def test_qcatch_remove_doublets() -> None:
    quant = await fixture_dir(TEST_DATA, "test_data/1k_pbmc_v3")
    result = await qcatch(quant, prefix="1k_pbmc_v3_doublets", chemistry="10X_3p_v3", remove_doublets=True)
    for f in (result.report, result.filtered_h5ad, result.metrics_summary):
        await assert_nonempty(f, label=f.path)


tests = [test_qcatch, test_qcatch_remove_doublets]
