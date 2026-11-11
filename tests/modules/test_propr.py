"""Tests for flyte_bio.modules.propr.

The cases and expected md5s are upstream's propr/propd and propr/grea tests
(modules commit 9a1cc9f), on the first contrast of upstream's contrast sheet
(``treatment_mCherry_hND6_``: treatment, mCherry versus hND6). As upstream,
the plot, session info, RDS and permutation FDRs are checked by name, and
grea's results (permutation-based) by name too.
"""

from flyte_bio.modules.propr import PropdResult, propr_grea, propr_propd
from tests.framework import assert_md5, env, fixture

EXPRESSION = "genomics/mus_musculus/rnaseq_expression/"
PREFIX = "treatment_mCherry_hND6_"
PROPD_VERSIONS_MD5 = "a6e655144d7d1d743781538442685666"
GENEWISE_MD5 = "3673576dbafb9e9f48ca707871f3f6a6"


async def run_propd(args: str) -> PropdResult:
    return await propr_propd(
        await fixture(EXPRESSION + "SRP254919.samplesheet.csv"),
        await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv"),
        "treatment",
        "mCherry",
        "hND6",
        prefix=PREFIX,
        args=args,
    )


async def check_common(result: PropdResult, genewise: str) -> None:
    await assert_md5(result.results_genewise, genewise, label="genewise")
    assert result.genewise_plot.name == f"{PREFIX}.propd.genewise.png"
    assert result.session_info.name == f"{PREFIX}.R_sessionInfo.log"
    await assert_md5(result.versions, PROPD_VERSIONS_MD5, label="versions.yml")


@env.task
async def test_propr_propd_defaults() -> None:
    # upstream case: "Test propr/propd using default params"
    result = await run_propd("--round_digits 5")
    await check_common(result, GENEWISE_MD5)
    optional = (result.rdata, result.results_pairwise, result.results_pairwise_filtered, result.adjacency, result.fdr)
    assert optional == (None,) * 5, "optional outputs written by default"


@env.task
async def test_propr_propd_all_outputs() -> None:
    # upstream case: "Test propr/propd when saving all outputs"
    result = await run_propd(
        "--save_pairwise_full true --save_pairwise true --save_adjacency true --save_rdata true --round_digits 5"
    )
    await check_common(result, GENEWISE_MD5)
    assert result.rdata is not None and result.rdata.name == f"{PREFIX}.propd.rds"
    assert result.results_pairwise and result.results_pairwise_filtered and result.adjacency
    await assert_md5(result.results_pairwise, "ff79479ea826c54aa178b0376784cbf1", label="pairwise")
    await assert_md5(result.results_pairwise_filtered, "9473f045ca1acfe1f9d6138528f743d1", label="pairwise_filtered")
    await assert_md5(result.adjacency, "3a947ffd8172990c6207aa7561a41cb8", label="adjacency")
    assert result.fdr is None, "FDRs written without permutations"


@env.task
async def test_propr_propd_boxcox() -> None:
    # upstream case: "Test propr/propd when using Box-cox transformation"
    result = await run_propd("--alpha 0.1 --round_digits 5")
    await check_common(result, "918f1bd489443235b2d0e22e12731985")


@env.task
async def test_propr_propd_permutation() -> None:
    # upstream case: "Test propr/propd when using permutation tests"
    result = await run_propd("--permutation 10 --number_of_cutoffs 10 --seed 123 --round_digits 5")
    await check_common(result, GENEWISE_MD5)
    assert result.fdr is not None and result.fdr.name == f"{PREFIX}.propd.fdr.tsv"


@env.task
async def test_propr_grea_chained_to_propd() -> None:
    # upstream case: "test grea chained to propd"
    propd = await run_propd("--round_digits 5 --save_adjacency true --features_id_col gene_name")
    assert propd.adjacency is not None, "propd wrote no adjacency matrix"
    gmt = await fixture("genomics/mus_musculus/gene_set_analysis/mh.all.v2022.1.Mm.symbols.gmt")
    result = await propr_grea(
        propd.adjacency, gmt, prefix=PREFIX, args="--permutation 10 --set_min 10 --seed 123 --round_digits 5"
    )
    assert result.results.name == f"{PREFIX}.grea.tsv"
    assert result.session_info.name == f"{PREFIX}.R_sessionInfo.log"
    await assert_md5(result.versions, "060fcd8ce4afc482e237fa75686a0aba", label="versions.yml")


tests = [
    test_propr_propd_defaults,
    test_propr_propd_all_outputs,
    test_propr_propd_boxcox,
    test_propr_propd_permutation,
    test_propr_grea_chained_to_propd,
]
