"""Tests for flyte_bio.modules.bff.

The cases and expected values are upstream's bff tests (modules commit
9a1cc9f): the HTO matrix with the combined and raw callers, with and without
preprocessing, and with ``--doHeatmap false`` (upstream's ``ext.args`` case).
As upstream, the matrix is staged as ``hto``, which the parameters table
records.
"""

from flyte.io import Dir

from flyte_bio.modules.bff import BffResult, bff
from tests.framework import assert_md5, env, fixture_dir

HTO = "genomics/homo_sapiens/10xgenomics/cellranger/hashing_demultiplexing/hto.tar.gz"
VERSIONS_MD5 = "71a50af83aebf56493d64a87f468a6f2"


async def bff_hto_matrix() -> Dir:
    return await fixture_dir(HTO, "hto")


async def check(result: BffResult, assignment: str, metrics: str, params: str) -> None:
    await assert_md5(result.assignment, assignment, label="assignment")
    await assert_md5(result.metrics, metrics, label="metrics")
    await assert_md5(result.params, params, label="params")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


@env.task
async def test_bff_combined_preprocessing() -> None:
    # upstream case: bff "hto_matrix - combined - preprocessing"
    result = await bff(await bff_hto_matrix(), prefix="test", methods="COMBINED", preprocessing=True)
    await check(
        result,
        "5346a38c04d6c6e7798f1b5e46c9ca91",
        "0b565a76802e8acf998bb35feeab35f8",
        "a8885b3d8aaa64b3350137e91856ca3b",
    )


@env.task
async def test_bff_combined_preprocessing_args() -> None:
    # upstream case: bff "hto_matrix - combined - preprocessing - test ext.args"
    result = await bff(
        await bff_hto_matrix(), prefix="test", methods="COMBINED", preprocessing=True, args="--doHeatmap false"
    )
    await check(
        result,
        "5346a38c04d6c6e7798f1b5e46c9ca91",
        "0b565a76802e8acf998bb35feeab35f8",
        "fbc8cbbe0c24be222cb2a0ea6261db99",
    )


@env.task
async def test_bff_raw_preprocessing() -> None:
    # upstream case: bff "hto_matrix - raw - preprocessing"
    result = await bff(await bff_hto_matrix(), prefix="test", methods="RAW", preprocessing=True)
    await check(
        result,
        "b9225b970358a61d62d2e67c0011b432",
        "14354e1e722db19e4dd9bbaabcce7ca5",
        "1e3c54b1dc9786e9e21c7baf2247ddef",
    )


@env.task
async def test_bff_combined_no_preprocessing() -> None:
    # upstream case: bff "hto_matrix - combined - no preprocessing"
    result = await bff(await bff_hto_matrix(), prefix="test", methods="COMBINED", preprocessing=False)
    await check(
        result,
        "d57dd3182d02b980ee231adfcef0581f",
        "0b565a76802e8acf998bb35feeab35f8",
        "a8885b3d8aaa64b3350137e91856ca3b",
    )


tests = [
    test_bff_combined_preprocessing,
    test_bff_combined_preprocessing_args,
    test_bff_raw_preprocessing,
    test_bff_combined_no_preprocessing,
]
