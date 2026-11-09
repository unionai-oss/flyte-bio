"""Tests for flyte_bio.modules.sparsesignatures.

The cases are upstream's sparsesignatures tests (modules commit 9a1cc9f,
SparseSignatures 2.16.0), with their ``ext.args`` (a small grid, one process,
seed 42). As upstream, the count matrix, the cross-validation MSEs and the
best configuration are checked by md5, the fit and plots for presence only,
and ``versions.yml`` by content (upstream snapshots its parsed YAML).
"""

from flyte.io import File

from flyte_bio.modules.sparsesignatures import SparsesignaturesResult, sparsesignatures
from tests.framework import assert_md5, assert_nonempty, env, fixture

ARGS = (
    "[K:2:3, nmf_runs:3, iterations:5, max_iterations_lasso:80, cross_validation_repetitions:3, "
    "cross_validation_iterations:3, cross_validation_entries:0.01, lambda_values_alpha:c(0.00, 0.01), "
    "lambda_values_beta:c(0.01, 0.05), lambda_rate_alpha:0, num_processes:NA, seed:42]"
)
VERSIONS = """\
"SPARSE_SIGNATURES":
    SparseSignatures: 2.16.0
    dplyr: 1.1.4
    ggplot2: 3.5.2
    patchwork: 1.3.0
    BSgenome.Hsapiens.1000genomes.hs37d5: 0.99.1
    BSgenome.Hsapiens.UCSC.hg38: 1.4.5
"""


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: SparsesignaturesResult, mut_counts: str, cv_means_mse: str, best_params: str) -> None:
    await assert_md5(result.mut_counts_rds, mut_counts, label=result.mut_counts_rds.name)
    await assert_md5(result.cv_means_mse_rds, cv_means_mse, label=result.cv_means_mse_rds.name)
    await assert_md5(result.best_params_config_rds, best_params, label=result.best_params_config_rds.name)
    for file in (result.nmf_lasso_out_rds, result.plot_pdf, result.plot_rds):
        await assert_nonempty(file, label=file.name)
    versions = await text(result.versions)
    assert versions == VERSIONS, f"versions.yml:\n{versions}\nexpected:\n{VERSIONS}"


@env.task
async def test_sparsesignatures_grch37() -> None:
    # upstream case: "test_sparsesignatures_GRCh37"
    table = await fixture("delete_me/sparsesignatures/test_GRCh37/MSeq_Set06_joint_table.tsv")
    result = await sparsesignatures([table], "GRCh37", prefix="test_grch37", args=ARGS)
    await check(
        result,
        "b4f87c7f3aca0211ff1371c03cdec4f8",
        "7f447130d4590f3f928857e5d48c0560",
        "fbf35a30b089d8a989815aada0dbcc94",
    )


@env.task
async def test_sparsesignatures_grch38() -> None:
    # upstream case: "test_sparsesignatures_GRCh38"
    table = await fixture("delete_me/sparsesignatures/test_GRCh38/HCC1395_WXS_HCC1395_joint_table.tsv")
    result = await sparsesignatures([table], "GRCh38", prefix="test_grch38", args=ARGS)
    await check(
        result,
        "78e5cdbfe5dbb34d6e1702b094ec05cf",
        "c9a8ebdd4f8bb61c5c0c366d412b1d6f",
        "39a02be6f840aea01de42e1d0cfba763",
    )


tests = [test_sparsesignatures_grch37, test_sparsesignatures_grch38]
