"""Tests for flyte_bio.modules.custom.

The cases and expected md5s are upstream's custom/filterdifferentialtable,
custom/matrixfilter, custom/tabulartogseagct, custom/tabulartogseacls and
custom/tabulartogseachip tests (modules commit 9a1cc9f), with their test
configs' ``args``.
"""

from flyte.io import File

from flyte_bio.modules.custom import (
    custom_filterdifferentialtable,
    custom_matrixfilter,
    custom_tabulartogseachip,
    custom_tabulartogseacls,
    custom_tabulartogseagct,
)
from tests.framework import assert_md5, env, fixture

EXPRESSION = "genomics/mus_musculus/rnaseq_expression/"


async def line_count(file: File) -> int:
    async with file.open("rb") as fh:
        return len(bytes(await fh.read()).decode().splitlines())


async def check_filtered(results, filtered: str, up: str, down: str, stat_column: str) -> None:
    await assert_md5(results.filtered, filtered, label="filtered")
    await assert_md5(results.filtered_up, up, label="filtered_up")
    await assert_md5(results.filtered_down, down, label="filtered_down")
    await assert_md5(results.versions, "c0d19fbbdd52d186345dd716bb94bd10", label="versions.yml")
    async with results.filtered.open("rb") as fh:
        header = bytes(await fh.read()).decode().splitlines()[0]
    assert "log2FoldChange" in header and stat_column in header, f"unexpected header {header!r}"
    # As upstream: the up and down tables partition the filtered rows (each with its own header).
    total = await line_count(results.filtered_up) + await line_count(results.filtered_down)
    assert total == await line_count(results.filtered) + 1, "up + down rows don't add up to the filtered rows"


@env.task
async def test_custom_filterdifferentialtable_padj() -> None:
    # upstream case: "Should filter differential expression table"
    table = await fixture(EXPRESSION + "SRP254919.salmon.merged.deseq2.results.tsv")
    results = await custom_filterdifferentialtable(table, "test", "log2FoldChange", 2, ">=", "padj", 0.05, "<=")
    await check_filtered(
        results,
        "f004047e216080aca0a4fc3ccde71179",
        "f004047e216080aca0a4fc3ccde71179",
        "24b3287c1e1162d047ced9dda147d54f",
        "padj",
    )


@env.task
async def test_custom_filterdifferentialtable_pvalue() -> None:
    # upstream case: "#lines(Up) + #lines(Down) = #lines(Full) + 1"
    table = await fixture(EXPRESSION + "SRP254919.salmon.merged.deseq2.results.tsv")
    results = await custom_filterdifferentialtable(table, "test", "log2FoldChange", 2, ">=", "pvalue", 0.05, "<=")
    await check_filtered(
        results,
        "75ff3b73def35be35592fe1ea8bf3bb4",
        "5b2db79939cab2f3a073f15a755d69b0",
        "ab76b7e861c6b07688825c1eb55d2da0",
        "pvalue",
    )


def matrixfilter_args(minimum_abundance: str) -> str:
    # As upstream's test params (the empty options leave double spaces, as upstream's join does).
    return (
        f'--sample_id_col "sample" --minimum_samples 1.0 --minimum_abundance {minimum_abundance}'
        '   --minimum_proportion_not_na "0.5"'
    )


async def run_matrixfilter(minimum_abundance: str, filtered: str, tests: str) -> None:
    abundance = await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.assay.tsv")
    samplesheet = await fixture(EXPRESSION + "SRP254919.samplesheet.sample_metadata.tsv")
    result = await custom_matrixfilter(abundance, "test", samplesheet, args=matrixfilter_args(minimum_abundance))
    await assert_md5(result.filtered, filtered, label="filtered")
    await assert_md5(result.tests, tests, label="tests")
    await assert_md5(result.versions, "275e382be9ca5bdd30ff5d6fc4340616", label="versions.yml")


@env.task
async def test_custom_matrixfilter_no_abundance_filter() -> None:
    # upstream case: "test_custom_matrixfilter - filtering_min_abundance false"
    await run_matrixfilter("false", "0c2224f6c33af1fe93c3b1382b3e05ce", "57f09670c8071ef449f3a78dab271681")


@env.task
async def test_custom_matrixfilter_min_abundance() -> None:
    # upstream case: "test_custom_matrixfilter"
    await run_matrixfilter("10", "ecc3ca4c6538f9248dce844108474a66", "bac204a5ac5f6bd2dcbd87297e16b686")


@env.task
async def test_custom_tabulartogseagct() -> None:
    # upstream case: "test_custom_tabulartogseagct"
    tabular = await fixture(EXPRESSION + "all.normalised_counts.tsv")
    gct = await custom_tabulartogseagct(tabular, "treatment_mCherry_hND6_")
    await assert_md5(gct, "9727aca45da7ebe619a1f964c9813eee", label="gct")


@env.task
async def test_custom_tabulartogseacls() -> None:
    # upstream case: "test_custom_tabulartogseacls"
    samples = await fixture(EXPRESSION + "SRP254919.samplesheet.sample_metadata.tsv")
    cls = await custom_tabulartogseacls(samples, "treatment_mCherry_hND6_sample_number", "treatment", separator="\t")
    await assert_md5(cls, "7b887159417842b1fe87afebe4f8c2b7", label="cls")


@env.task
async def test_custom_tabulartogseachip() -> None:
    # upstream case: "test tsv to chip"
    tabular = await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv")
    chip = await custom_tabulartogseachip(
        tabular, "SRP254919.salmon.merged.gene_counts.top1000cov", "gene_id", "gene_name"
    )
    await assert_md5(chip, "2ab8a685c675ce2fb97142526766044a", label="chip")


tests = [
    test_custom_filterdifferentialtable_padj,
    test_custom_filterdifferentialtable_pvalue,
    test_custom_matrixfilter_no_abundance_filter,
    test_custom_matrixfilter_min_abundance,
    test_custom_tabulartogseagct,
    test_custom_tabulartogseacls,
    test_custom_tabulartogseachip,
]
