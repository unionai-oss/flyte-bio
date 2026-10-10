"""Tests for flyte_bio.modules.qcatch.

Upstream's test data: a simpleaf quantification of 1k PBMCs (10x v3) from
QCatch's authors. Upstream snapshots the md5 of the filtered h5ad and of the
metrics summary; the doublet run only checks that QCatch succeeds.
"""

from flyte_bio.modules.qcatch import qcatch
from tests.framework import assert_md5, assert_nonempty, env, fixture_dir

TEST_DATA = "https://umd.box.com/shared/static/zd4sai70uw9fs24e1qx6r41ec50pf45g.zip"


@env.task
async def test_qcatch() -> None:
    quant = await fixture_dir(TEST_DATA, "test_data/1k_pbmc_v3")
    result = await qcatch(quant, prefix="1k_pbmc_v3", chemistry="10X_3p_v3")
    assert result.report.path.endswith("/1k_pbmc_v3_qcatch_report.html"), result.report.path
    await assert_nonempty(result.report, label="report")
    # The h5ad's bytes vary between runs (HDF5), so unlike upstream only its existence is checked;
    # the metrics summary, which holds the cell count and QC metrics, is pinned instead.
    await assert_nonempty(result.filtered_h5ad, label="filtered h5ad")
    await assert_md5(result.metrics_summary, "22deb2f3aea4e4abf476970576362428", label="metrics summary")


@env.task
async def test_qcatch_remove_doublets() -> None:
    quant = await fixture_dir(TEST_DATA, "test_data/1k_pbmc_v3")
    result = await qcatch(quant, prefix="1k_pbmc_v3_doublets", chemistry="10X_3p_v3", remove_doublets=True)
    for f in (result.report, result.filtered_h5ad, result.metrics_summary):
        await assert_nonempty(f, label=f.path)


tests = [test_qcatch, test_qcatch_remove_doublets]
