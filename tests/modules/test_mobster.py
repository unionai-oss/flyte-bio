"""Tests for flyte_bio.modules.mobster.

The cases are upstream's mobster tests (modules commit 9a1cc9f), with their
test config's ``ext.args`` (the module's default). The fits start from random
initial conditions, so, as upstream, the RDS, PDF and PNG outputs are checked
for presence only, and ``versions.yml`` by md5.
"""

from flyte_bio.modules.mobster import MobsterResult, mobster
from tests.framework import assert_md5, assert_nonempty, env, fixture

DATA = "delete_me/mobster/HCC1395_WXS_HCC1395_"
VERSIONS_MD5 = "d7bc42ac0521d8f55d715eab700aa3f4"


async def check(result: MobsterResult) -> None:
    for file in (
        result.fit_rds,
        result.best_fit_rds,
        result.best_plots_rds,
        result.report_rds,
        result.report_pdf,
        result.report_png,
    ):
        await assert_nonempty(file, label=file.name)
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


@env.task
async def test_mobster_rds_input() -> None:
    # upstream case: "test rds input"
    joint = await fixture(DATA + "multi_cnaqc_ALL.rds")
    await check(await mobster(joint, patient="P1", tumour_sample="HCC1395_HCC1395T", prefix="test"))


@env.task
async def test_mobster_csv_input() -> None:
    # upstream case: "test csv input"
    joint = await fixture(DATA + "joint_table.csv")
    await check(await mobster(joint, patient="P1", tumour_sample="HCC1395_HCC1395T", prefix="test"))


tests = [test_mobster_rds_input, test_mobster_csv_input]
