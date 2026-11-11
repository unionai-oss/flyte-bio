"""Tests for flyte_bio.modules.gprofiler2.

The cases and expected md5s are upstream's gprofiler2/gost tests (modules
commit 9a1cc9f), with their test params' ``args`` (organism mode, so
g:Profiler's web service is queried, so the md5s hold while g:Profiler serves
the same data release). As upstream, the plots, HTML, RDS and session info
are checked by name.
"""

from flyte.io import File

from flyte_bio.modules.gprofiler2 import GostResult, gprofiler2_gost
from tests.framework import assert_md5, env, fixture

EXPRESSION = "genomics/mus_musculus/rnaseq_expression/"
PREFIX = "Condition_genotype_WT_KO"
BASE = PREFIX + ".gprofiler2."


def gost_args(sources: str | None) -> str:
    # As upstream's test params, joined with their empty (null) options.
    options = [
        '--significant "true"',
        '--measure_underrepresentation "false"',
        '--correction_method "gSCS"',
        '--evcodes "false"',
        '--pval_threshold "0.05"',
        "--domain_scope annotated",
        '--min_diff "1"',
        "--round_digits 4",
        '--palette_name "Blues"',
        '--de_id_column "gene_id"',
        "",
        '--organism "mmusculus"',
        "",
        "" if sources is None else f'--sources "{sources}"',
    ]
    return " ".join(options).strip()


async def run(sources: str | None) -> GostResult:
    return await gprofiler2_gost(
        await fixture(EXPRESSION + "Condition_genotype_WT_KO.deseq2.results_filtered.tsv"),
        prefix=PREFIX,
        gmt_file=await fixture("genomics/mus_musculus/gene_set_analysis/mh.all.v2022.1.Mm.symbols.gmt"),
        background_file=await fixture(EXPRESSION + "study.filtered.tsv"),
        args=gost_args(sources),
    )


async def check_common(result: GostResult, sub: dict[str, str], all_md5: str) -> None:
    await assert_md5(result.all_enrich, all_md5, label="all_enriched_pathways")
    names = [f"{BASE}{source}.sub_enriched_pathways.tsv" for source in sub]
    assert [f.name for f in result.sub_enrich] == names, [f.name for f in result.sub_enrich]
    for file, md5 in zip(result.sub_enrich, sub.values(), strict=True):
        await assert_md5(file, md5, label=file.name)
    plots = [f"{BASE}{source}.sub_enriched_pathways.png" for source in sub]
    assert [f.name for f in result.sub_plot] == plots, [f.name for f in result.sub_plot]
    assert result.plot_png is not None and result.plot_png.name == BASE + "gostplot.png"
    assert result.plot_html is not None and result.plot_html.name == BASE + "gostplot.html"
    assert result.rds is not None and result.rds.name == BASE + "gost_results.rds"
    assert result.session_info.name == "R_sessionInfo.log"
    async with result.versions.open("rb") as fh:
        versions = bytes(await fh.read()).decode()
    # Upstream's snapshot records g:Profiler data build e114_eg62_p19_2fe5c4c1; the service has since moved to
    # another build of the same release (e114_eg62_p19_9dc92ac8), with the same results, so only the release is
    # checked.
    for line in ("r-base: 4.3.3", "r-ggplot2: 3.4.3", "r-gprofiler2: 0.2.2", "gprofiler-data: e114_eg62_p19_"):
        assert line in versions, f"versions.yml lacks {line!r}: {versions}"


@env.task
async def test_gprofiler2_gost() -> None:
    # upstream case: "test"
    result = await run(None)
    sub = {
        "GO:BP": "1201012160e769c4e1f0d6fd1a20e618",
        "GO:CC": "06d30deaee5c3ba5f7bd0922d683bb3d",
        "GO:MF": "036e130863cc281d3771e5f114119c42",
        "HP": "6ded662a2f00ca35b09df9b81cf89ade",
        "KEGG": "dd88cd34b1db98eb47acb123a85e68f2",
        "REAC": "6328990ca27a958621d1d468b880c1a8",
        "TF": "d3ad0c6a9a9cd01f2e9cf0bbf2e2fa69",
        "WP": "bf45e6db392a678b7fd5c6b277fbd58d",
    }
    await check_common(result, sub, "17c840d9a837ec0e3f60ca519f0956ee")
    assert result.filtered_gmt is None, "a filtered GMT without --sources"


@env.task
async def test_gprofiler2_gost_sources() -> None:
    # upstream case: "test-sources"
    result = await run("GO:CC,WP")
    sub = {"GO:CC": "06d30deaee5c3ba5f7bd0922d683bb3d", "WP": "bf45e6db392a678b7fd5c6b277fbd58d"}
    await check_common(result, sub, "4eca2a382d9e5452cb54da925a66357a")
    gmt: File | None = result.filtered_gmt
    assert gmt is not None and gmt.name == "gprofiler_full_mmusculus.GO:CC_WP.ENSG_filtered.gmt"
    await assert_md5(gmt, "fe7fb3b351e0ffa5567f866d103c5331", label="filtered GMT")


tests = [test_gprofiler2_gost, test_gprofiler2_gost_sources]
