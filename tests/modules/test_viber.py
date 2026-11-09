"""Tests for flyte_bio.modules.viber.

The cases are upstream's viber tests (modules commit 9a1cc9f), with their test
config's ``ext.args`` (the module's default). VIBER's fit is unseeded, so, as
upstream, the RDS, PDF and PNG outputs are checked for presence only, and
``versions.yml`` by md5.
"""

from flyte.io import File

from flyte_bio.modules.viber import ViberResult, viber
from tests.framework import assert_md5, assert_nonempty, env, fixture

VERSIONS_MD5 = "0aff4c2ae3de260cbef5870e5f4f0374"


async def check(result: ViberResult) -> None:
    for file in (
        result.best_fit_rds,
        result.best_heuristic_fit_rds,
        result.best_fit_plots_rds,
        result.best_heuristic_fit_plots_rds,
        result.report_rds,
        result.report_png,
        result.report_pdf,
    ):
        await assert_nonempty(file, label=file.name)
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


async def run_viber(joint: File) -> ViberResult:
    return await viber(joint, patient="HCC1395", tumour_samples=["HCC1395_HCC1395T"], prefix="test")


@env.task
async def test_viber_rds() -> None:
    # upstream case: "viber - rds"
    await check(await run_viber(await fixture("delete_me/viber/HCC1395_WXS_HCC1395_multi_cnaqc_ALL.rds")))


@env.task
async def test_viber_csv() -> None:
    # upstream case: "viber - csv"
    await check(await run_viber(await fixture("delete_me/mobster/HCC1395_WXS_HCC1395_joint_table.csv")))


tests = [test_viber_rds, test_viber_csv]
