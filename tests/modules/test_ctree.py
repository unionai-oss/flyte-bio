"""Tests for flyte_bio.modules.ctree.

The first two cases are upstream's ctree tests (modules commit 9a1cc9f, ctree
1.1.0): a MOBSTER best fit and a PyClone-VI cluster table, checked as
upstream by the outputs' names and ``versions.yml``'s md5. A MOBSTER fit's
outputs go under its sample's subdirectory.

The chain cases feed ctree the outputs of this package's own MOBSTER and
PyClone-VI tasks, run on those modules' test inputs, as the tumour evolution
pipeline does.
"""

from flyte_bio.modules.ctree import CtreeResult, ctree
from flyte_bio.modules.mobster import mobster
from flyte_bio.modules.pyclonevi import pyclonevi
from tests.framework import assert_md5, assert_nonempty, env, fixture

VERSIONS_MD5 = "fd370aeddbd6a0c86437152d14919289"


async def check(result: CtreeResult, tool: str, subdir: str = "") -> None:
    expected = {
        "trees_rds": f"test_ctree_{tool}.rds",
        "plots_rds": f"test_ctree_{tool}_plots.rds",
        "report_rds": f"test_ctree_{tool}_report.rds",
        "report_pdf": f"test_ctree_{tool}_report.pdf",
        "report_png": f"test_ctree_{tool}_report.png",
    }
    for field, name in expected.items():
        file = getattr(result, field)
        assert file is not None, f"ctree wrote no {name}"
        assert file.name == name, f"ctree {field}: {file.name}, expected {name}"
        await assert_nonempty(file, label=name)
        if subdir:
            assert await result.results.get_file(f"{subdir}/{name}") is not None, f"{name} isn't under {subdir}/"
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


@env.task
async def test_ctree_mobster_rds_input() -> None:
    # upstream case: "test mobster rds input"
    fit = await fixture("delete_me/ctree/HCC1395_WXS_HCC1395_HCC1395_HCC1395T_mobster_st_best_fit.rds")
    await check(await ctree(fit, prefix="test"), "mobster", subdir="HCC1395_HCC1395T")


@env.task
async def test_ctree_pyclonevi_tsv_input() -> None:
    # upstream case: "test pyclonevi tsv input"
    table = await fixture("delete_me/ctree/HCC1395_WXS_HCC1395_cluster_table.csv")
    await check(await ctree(table, prefix="test"), "pyclonevi")


@env.task
async def test_ctree_mobster_chain() -> None:
    joint = await fixture("delete_me/mobster/HCC1395_WXS_HCC1395_joint_table.csv")
    fit = await mobster(joint, patient="P1", tumour_sample="HCC1395_HCC1395T", prefix="test")
    await check(await ctree(fit.best_fit_rds, prefix="test"), "mobster", subdir="HCC1395_HCC1395T")


@env.task
async def test_ctree_pyclonevi_chain() -> None:
    joint = await fixture("delete_me/pyclonevi/HCC1395_WXS_HCC1395_joint_table.tsv")
    clusters = await pyclonevi(joint, patient="HCC1395", tumour_samples=["HCC1395_HCC1395T"], prefix="test")
    await check(await ctree(clusters.cluster_table, prefix="test"), "pyclonevi")


tests = [
    test_ctree_mobster_rds_input,
    test_ctree_pyclonevi_tsv_input,
    test_ctree_mobster_chain,
    test_ctree_pyclonevi_chain,
]
