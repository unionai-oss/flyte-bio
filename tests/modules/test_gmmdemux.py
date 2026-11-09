"""Tests for flyte_bio.modules.gmmdemux.

The cases are upstream's gmmdemux tests (modules commit 9a1cc9f): the HTO
matrix with hashtags ``MS-11,MS-12``, with the simplified report, and with
the full and summary reports and the test config's ``ext.args``. As upstream,
the outputs are checked for presence only.
"""

from flyte_bio.modules.gmmdemux import GmmDemuxResult, gmmdemux
from tests.framework import assert_nonempty, env, fixture_dir

HTO = "genomics/homo_sapiens/10xgenomics/cellranger/hashing_demultiplexing/hto.tar.gz"


async def check(result: GmmDemuxResult) -> None:
    for f in (result.matrix, result.barcodes, result.features, result.config_report):
        await assert_nonempty(f, label=f.name)


@env.task
async def test_gmmdemux_simple() -> None:
    # upstream case: gmmdemux "Standard_Multiome - 10x mtx - simple"
    result = await gmmdemux(await fixture_dir(HTO, "hto"), "MS-11,MS-12", prefix="test")
    await check(result)


@env.task
async def test_gmmdemux_full_summary() -> None:
    # upstream case: gmmdemux "Standard_Multiome - 10x mtx"
    result = await gmmdemux(
        await fixture_dir(HTO, "hto"),
        "MS-11,MS-12",
        prefix="test",
        full_report=True,
        summary_report=True,
        args="-t 0.8 --summary 11089 -a 0.05",
    )
    await check(result)
    assert result.classification_report.name == "GMM_full.csv", result.classification_report.name
    assert result.summary_report is not None, "no summary report"


tests = [test_gmmdemux_simple, test_gmmdemux_full_summary]
