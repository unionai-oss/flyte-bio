"""Tests for flyte_bio.modules.scanpy's HashSolo.

The cases and expected values are upstream's scanpy/hashsolo tests (modules
commit 9a1cc9f): the HTO matrix (staged as ``hto``) and upstream's generated
h5ad with hashing columns ``0`` to ``9``, each with ``--round_digits 10`` and
with priors ``0.20 0.8 0.19`` too (upstream's test config ``ext.args``). As
upstream, the parameters table and versions file are compared by md5 and the
assignment and h5ad (not stable) checked by name.
"""

from flyte.io import File

from flyte_bio.modules.scanpy import HashsoloResult, scanpy_hashsolo
from tests.framework import assert_md5, assert_nonempty, env, fixture, fixture_dir

DATA = "genomics/homo_sapiens/10xgenomics/cellranger/hashing_demultiplexing/"
VERSIONS_MD5 = "f6178c11c8b15300ccafa689830c64e3"
ROUND = "--round_digits 10"
PRIORS = "--priors 0.20 0.8 0.19 --round_digits 10"
COLUMNS = [str(i) for i in range(10)]


async def scanpy_hashsolo_h5ad() -> File:
    return await fixture(DATA + "hashsolo_anndata.h5ad")


async def check(result: HashsoloResult, params: str) -> None:
    await assert_md5(result.params, params, label="params")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")
    for f, name in ((result.assignment, "test_assignment_hashsolo.csv"), (result.h5ad, "test_hashsolo.h5ad")):
        assert f.name == name, f"expected {name}, got {f.name}"
        await assert_nonempty(f, label=name)


@env.task
async def test_scanpy_hashsolo_10x() -> None:
    # upstream case: scanpy/hashsolo "hto_data_10x"
    result = await scanpy_hashsolo("test", matrix=await fixture_dir(DATA + "hto.tar.gz", "hto"), args=ROUND)
    await check(result, "17ed9d21353013b8547ac829874ff514")


@env.task
async def test_scanpy_hashsolo_10x_args() -> None:
    # upstream case: scanpy/hashsolo "hto_data_10x - ext.args"
    result = await scanpy_hashsolo("test", matrix=await fixture_dir(DATA + "hto.tar.gz", "hto"), args=PRIORS)
    await check(result, "4b992e5daf73b8f23f5d847df32dc177")


@env.task
async def test_scanpy_hashsolo_h5ad() -> None:
    # upstream case: scanpy/hashsolo "hto_data_h5ad (generated, no real data)"
    result = await scanpy_hashsolo("test", h5ad=await scanpy_hashsolo_h5ad(), cell_hashing_columns=COLUMNS, args=ROUND)
    await check(result, "7054d81f10cb650a72f9dea5bb23f44b")


@env.task
async def test_scanpy_hashsolo_h5ad_args() -> None:
    # upstream case: scanpy/hashsolo "hto_data_h5ad (generated, no real data) - ext.args"
    result = await scanpy_hashsolo("test", h5ad=await scanpy_hashsolo_h5ad(), cell_hashing_columns=COLUMNS, args=PRIORS)
    await check(result, "c9e859256e6d28dc03049496b7eee509")


tests = [
    test_scanpy_hashsolo_10x,
    test_scanpy_hashsolo_10x_args,
    test_scanpy_hashsolo_h5ad,
    test_scanpy_hashsolo_h5ad_args,
]
