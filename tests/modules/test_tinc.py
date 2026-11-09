"""Tests for flyte_bio.modules.tinc.

The case is upstream's tinc test (modules commit 9a1cc9f, TINC 0.1.0), with
its ``ext.args`` (the module's default). As upstream, the fit, plots and QC
table are checked for presence only, and ``versions.yml`` by md5.
"""

from flyte_bio.modules.tinc import tinc
from tests.framework import assert_md5, assert_nonempty, env, fixture

DATA = "delete_me/tinc/HCC1395_WXS_HCC1395_HCC1395_HCC1395T_"


@env.task
async def test_tinc_rds() -> None:
    # upstream case: "tinc - rds"
    cna = await fixture(DATA + "cna.rds")
    snv = await fixture(DATA + "snv.rds")
    result = await tinc(cna, snv, tumour_sample="HCC1395_HCC1395T", normal_sample="HCC1395_HCC1395N", prefix="test")
    for file in (result.fit_rds, result.plot_rds, result.plot_pdf, result.qc_csv):
        await assert_nonempty(file, label=file.name)
    await assert_md5(result.versions, "fd8574ef3613bda9da0db50875ab36fc", label="versions.yml")


tests = [test_tinc_rds]
