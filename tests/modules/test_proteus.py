"""Tests for flyte_bio.modules.proteus.

The case and expected md5s are upstream's proteus/readproteingroups test
(modules commit 9a1cc9f), with its test params' ``args``. As upstream, the
plots, RDS files and session info are checked by name.
"""

from flyte_bio.modules.proteus import proteus_readproteingroups
from tests.framework import assert_md5, env, fixture

MAXQUANT = "proteomics/maxquant/"


@env.task
async def test_proteus_readproteingroups() -> None:
    # upstream case: "test_proteus_readproteingroups"
    result = await proteus_readproteingroups(
        await fixture(MAXQUANT + "MaxQuant_samplesheet.tsv"),
        await fixture(MAXQUANT + "MaxQuant_proteinGroups.txt"),
        args='--contrast_variable fakeBatch --sample_id_col Experiment  --measure_col_prefix "LFQ intensity" '
        '--round_digits 3 --norm_function "normalizeQuantiles"',
    )
    names = {
        result.dendro_plot.name: "normalizeQuantiles.normalized_dendrogram.png",
        result.mean_var_plot.name: "normalizeQuantiles.normalized_mean_variance_relationship.png",
        result.raw_dist_plot.name: "raw_distributions.png",
        result.norm_dist_plot.name: "normalizeQuantiles.normalized_distributions.png",
        result.raw_rdata.name: "raw_proteingroups.rds",
        result.norm_rdata.name: "normalizeQuantiles.normalized_proteingroups.rds",
        result.session_info.name: "R_sessionInfo.log",
    }
    for got, expected in names.items():
        assert got == expected, f"{got}, expected {expected}"
    await assert_md5(result.raw_tab, "af7a556997b82f4ba38bfae8a6b91301", label="raw table")
    assert result.norm_tab.name == "normalizeQuantiles.normalized_proteingroups_tab.tsv"
    await assert_md5(result.norm_tab, "4a0a74febfc6ba921210f465e890ff5b", label="normalised table")
    await assert_md5(result.versions, "e05128fdfdbcf466f3d9307c6493e167", label="versions.yml")


tests = [test_proteus_readproteingroups]
