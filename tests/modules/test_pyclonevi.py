"""Tests for flyte_bio.modules.pyclonevi.

The case is upstream's pyclonevi test (modules commit 9a1cc9f, PyClone-VI
0.2.0), with its test config's ``ext.args`` (the module's default,
``[n_cluster:20]``). PyClone-VI's restarts are unseeded, so, as upstream, the
outputs are checked by name only (and here for being nonempty), and
``versions.yml`` by md5.
"""

from flyte_bio.modules.pyclonevi import pyclonevi
from tests.framework import assert_md5, assert_nonempty, env, fixture


@env.task
async def test_pyclonevi_run() -> None:
    # upstream case: "pyclonevi"
    joint = await fixture("delete_me/pyclonevi/HCC1395_WXS_HCC1395_joint_table.tsv")
    result = await pyclonevi(joint, patient="HCC1395", tumour_samples=["HCC1395_HCC1395T"], prefix="test")
    names = [result.cluster_table.name, *(f.name for f in result.pyclone_input)]
    names += [result.all_fits.name, result.best_fit.name]
    expected = [
        "test_cluster_table.csv",
        "test_pyclone_input.tsv",
        "test_pyclone_input_all_samples.tsv",
        "test_all_fits.h5",
        "test_best_fit.txt",
    ]
    assert names == expected, f"pyclonevi outputs {names}, expected {expected}"
    for file in (result.cluster_table, *result.pyclone_input, result.all_fits, result.best_fit):
        await assert_nonempty(file, label=file.name)
    await assert_md5(result.versions, "3d517410bfd758e0eea9bef7f40d22e6", label="versions.yml")


tests = [test_pyclonevi_run]
