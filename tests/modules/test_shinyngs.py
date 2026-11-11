"""Tests for flyte_bio.modules.shinyngs.

The cases and expected values are upstream's shinyngs/validatefomcomponents,
shinyngs/staticexploratory, shinyngs/staticdifferential and shinyngs/app
tests (modules commit 9a1cc9f), with their ``args``. As upstream, plots and
the app's ``data.rds`` are checked by name, the validated tables and
``app.R`` by md5. The app cases copy the count matrix to ``normalised.tsv``
and the results table to ``second_contrast_stats.tsv``, as upstream's do.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.shinyngs import (
    shinyngs_app,
    shinyngs_staticdifferential,
    shinyngs_staticexploratory,
    shinyngs_validatefomcomponents,
)
from tests.framework import assert_md5, env, fixture

EXPRESSION = "genomics/mus_musculus/rnaseq_expression/"
GENE_SET_ANALYSIS = "genomics/mus_musculus/gene_set_analysis/"
APP_R_MD5 = "1e02e7b6d263199c3f45891b3126a3ca"


async def copy_as(file: File, name: str) -> File:
    async with file.open("rb") as fh:
        data = bytes(await fh.read())
    path = Path(tempfile.mkdtemp()) / name
    path.write_bytes(data)
    return await File.from_local(str(path))


async def run_validate(contrasts: str, contrasts_md5: str) -> None:
    result = await shinyngs_validatefomcomponents(
        await fixture(EXPRESSION + "SRP254919.samplesheet.csv"),
        [await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv")],
        await fixture(EXPRESSION + contrasts),
        prefix="SRP254919",
        feature_meta=await fixture(EXPRESSION + "Mus_musculus.anno.tsv"),
        args="--sample_id_col 'sample' --feature_id_col 'gene_id'",
    )
    assert result.sample_meta.name == "SRP254919.samplesheet.sample_metadata.tsv"
    await assert_md5(result.sample_meta, "7c5233cf34fd8eaf7071f21f6fcac7d8", label="sample metadata")
    assert result.feature_meta is not None and result.feature_meta.name == "Mus_musculus.anno.feature_metadata.tsv"
    await assert_md5(result.feature_meta, "1b6a2af23f009f4669401678fe0947bb", label="feature metadata")
    assert [f.name for f in result.assays] == ["SRP254919.salmon.merged.gene_counts.top1000cov.assay.tsv"]
    await assert_md5(result.assays[0], "0c2224f6c33af1fe93c3b1382b3e05ce", label="assay")
    assert result.contrasts.name == "SRP254919.contrasts.contrasts_file.tsv"
    await assert_md5(result.contrasts, contrasts_md5, label="contrasts")


@env.task
async def test_shinyngs_validatefomcomponents() -> None:
    # upstream case: "test"
    await run_validate("SRP254919.contrasts.csv", "e0980240404f48216aedb63b5094029f")


@env.task
async def test_shinyngs_validatefomcomponents_yaml() -> None:
    # upstream case: "test_yaml"
    await run_validate("SRP254919.contrasts.yaml", "0329c0b0b811e3c218453ae70b25496e")


async def run_exploratory(args: str):
    return await shinyngs_staticexploratory(
        await fixture(EXPRESSION + "SRP254919.samplesheet.csv"),
        await fixture(EXPRESSION + "SRP254919.gene_meta.tsv"),
        [await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv")],
        "treatment",
        prefix="treatment",
        args=args,
    )


def check_png_names(result) -> None:
    pngs = (result.boxplots_png, result.densities_png, result.mad_png, result.pca2d_png, result.pca3d_png)
    names = [f.name if f else None for f in (*pngs, result.dendro)]
    expected = ["boxplot.png", "density.png", "mad_correlation.png", "pca2d.png", "pca3d.png", "sample_dendrogram.png"]
    assert names == expected, names


@env.task
async def test_shinyngs_staticexploratory() -> None:
    # upstream case: "mouse - defaults"
    check_png_names(await run_exploratory(""))


@env.task
async def test_shinyngs_staticexploratory_specify_log() -> None:
    # upstream case: "mouse - specify log"
    check_png_names(await run_exploratory("--log2_assays '1'"))


@env.task
async def test_shinyngs_staticexploratory_html() -> None:
    # upstream case: "mouse - html". Upstream's snapshot records that none of the HTML files holds its
    # probe lines as a whole line (readLines().contains is false for each); this checks the same.
    result = await run_exploratory("--write_html")
    check_png_names(result)
    probes = {
        "boxplots_html": 'ENSMUSG00000027456","Gm37080',
        "densities_html": "-1.4916353753463203,-1.3078530974016225,-1.1240708194569247,-0.94028854151222707",
        "mad_html": "0,-0.74295280067699376,0.67449075947659531,-0.6744907594765952",
        "pca2d_html": 'SRX8042381","SRX8042382',
        "pca3d_html": 'SRX8042381","SRX8042382',
    }
    for output, line in probes.items():
        file = getattr(result, output)
        assert file is not None, f"no {output}"
        async with file.open("rb") as fh:
            lines = bytes(await fh.read()).decode().splitlines()
        assert line not in lines, f"{output} unexpectedly holds the line {line!r}"


@env.task
async def test_shinyngs_staticdifferential() -> None:
    # upstream case: "test"
    args = (
        '--feature_id_col "gene_id" --reference_level "mCherry" --treatment_level "hND6" '
        '--fold_change_col "log2FoldChange" --p_value_column "padj" --diff_feature_id_col "gene_id" '
        '--fold_change_threshold "2.0" --p_value_threshold "0.05" --unlog_foldchanges --palette_name "Set1"'
    )
    result = await shinyngs_staticdifferential(
        await fixture(EXPRESSION + "treatment_mCherry_hND6_.deseq2.results.tsv"),
        await fixture(EXPRESSION + "Mus_musculus.anno.tsv"),
        prefix="treatment_mCherry_hND6_",
        args=args,
    )
    assert result.volcanos_png.name == "volcano.png"
    assert result.volcanos_html is None, "volcano.html without --write_html"


async def run_app(multi_matrix: bool, args: str, enrichment: bool = False) -> None:
    raw = await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv")
    differential = await fixture(EXPRESSION + "SRP254919.salmon.merged.deseq2.results.tsv")
    assays = [raw, await copy_as(raw, "normalised.tsv")] if multi_matrix else [raw]
    gene_sets, enrichment_results = [], []
    if enrichment:
        gene_sets = [await fixture(GENE_SET_ANALYSIS + "mh.all.v2022.1.Mm.symbols.gmt")]
        enrichment_results = [
            await fixture(
                f"{GENE_SET_ANALYSIS}shinyngs/{contrast}.mh.all.v2022.1.Mm.symbols.gsea_report_for_{level}.tsv"
            )
            for contrast in ("treatment_mCherry_hND6_", "treatment_mCherry_hND6_sample_number")
            for level in ("hND6", "mCherry")
        ]
    result = await shinyngs_app(
        await fixture(EXPRESSION + "SRP254919.samplesheet.csv"),
        await fixture(EXPRESSION + "SRP254919.gene_meta.tsv"),
        assays,
        await fixture(EXPRESSION + "SRP254919.contrasts.csv"),
        [differential, await copy_as(differential, "second_contrast_stats.tsv")],
        2 if multi_matrix else 1,
        prefix="SRP254919_test",
        gene_sets=gene_sets,
        enrichment_results=enrichment_results,
        args=args,
    )
    assert result.data.name == "data.rds"
    await assert_md5(result.app, APP_R_MD5, label="app.R")


@env.task
async def test_shinyngs_app_multi_matrix() -> None:
    # upstream case: "mouse - multi matrix"
    await run_app(True, "--assay_names raw,normalised")


@env.task
async def test_shinyngs_app_enrichment() -> None:
    # upstream case: "mouse - with enrichment"
    await run_app(
        True,
        "--assay_names raw,normalised --enrichment_gene_sets mh.all.v2022.1.Mm.symbols.gmt "
        "--enrichment_filename_template '{contrast_name}.{geneset_type}.gsea_report_for_{target|reference}.tsv' "
        "--enrichment_gene_type_id gene_name",
        enrichment=True,
    )


@env.task
async def test_shinyngs_app_single_matrix() -> None:
    # upstream case: "mouse - single matrix"
    await run_app(False, "")


tests = [
    test_shinyngs_validatefomcomponents,
    test_shinyngs_validatefomcomponents_yaml,
    test_shinyngs_staticexploratory,
    test_shinyngs_staticexploratory_specify_log,
    test_shinyngs_staticexploratory_html,
    test_shinyngs_staticdifferential,
    test_shinyngs_app_multi_matrix,
    test_shinyngs_app_enrichment,
    test_shinyngs_app_single_matrix,
]
