"""Tests for flyte_bio.modules.limma.

The cases and expected values are upstream's limma/differential tests
(modules commit 9a1cc9f), with each case's test config ``args``: the model
(and, where upstream snapshots them, results and normalised counts) by md5,
and the result values upstream checks for. The microarray cases normalise
GSE38751's CEL files with affy's justRMA first, as upstream's setup does.
"""

from flyte.io import File

from flyte_bio.modules.affy import affy_justrma
from flyte_bio.modules.limma import LimmaDifferentialResult, limma_differential
from tests.framework import assert_md5, env, fixture, fixture_dir

DREAM = "https://github.com/nf-core/test-datasets/raw/differentialabundance/modules_testdata/variancepartition_dream/"
ARRAY = "genomics/homo_sapiens/array_expression/"
EXPRESSION = "genomics/mus_musculus/rnaseq_expression/"
VERSIONS_MD5 = "88a6e42d753077edab8daf829cd4d943"
ARRAY_MODEL = "70b000f632b8bdba4917046362dd876b"
SHEET_MODEL = "d2113d82b76046c319e6602da2ad74d6"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: LimmaDifferentialResult, model: str, contains: dict[str, list[str]] | None = None) -> None:
    await assert_md5(result.model, model, label="model")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")
    for output, values in (contains or {}).items():
        file = getattr(result, output)
        assert file is not None, f"no {output}"
        body = await text(file)
        for value in values:
            assert value in body, f"{output} lacks {value!r}"


@env.task
async def test_limma_formula_interaction() -> None:
    # upstream case: "RNAseq - Feature Counts - formula + comparison contrast string - interaction"
    result = await limma_differential(
        await fixture(DREAM + "metadata.tsv"),
        await fixture(DREAM + "counts.tsv"),
        prefix="genotype_WT_KO_treatment_Control_Treated",
        formula="~ genotype * treatment",
        comparison="genotypeWT.treatmentTreated",
        args="--sample_id_col sample --probe_id_col gene_id",
    )
    await check(result, "82165b8b24da8f8cffafa457290feeb3")


async def run_array(args: str, **contrast: str) -> LimmaDifferentialResult:
    """limma on upstream's justRMA matrix of GSE38751 (contrast ID ``diagnosis_normal_uremia``)."""
    sheet = await fixture(ARRAY + "GSE38751.csv")
    rma = await affy_justrma(
        sheet, await fixture_dir(ARRAY + "GSE38751_RAW.tar"), "test", args="--sample_name_col name"
    )
    return await limma_differential(sheet, rma.expression, prefix="diagnosis_normal_uremia", args=args, **contrast)


DIAGNOSIS = {"contrast_variable": "diagnosis", "reference": "normal", "target": "uremia"}


@env.task
async def test_limma_array() -> None:
    # upstream case: "test_limma_differential"
    result = await run_array("--sample_id_col name --blocking_variables null", **DIAGNOSIS)
    await check(
        result,
        ARRAY_MODEL,
        {"session_info": ["limma_3.58.1"], "results": ["1007_s_at\t-0.2775254", "1053_at\t-0.071547786"]},
    )


@env.task
async def test_limma_array_formula() -> None:
    # upstream case: "test_limma_differential - with formula and with complex contrasts"
    result = await run_array(
        "--sample_id_col name --blocking_variables null", formula="~ diagnosis", comparison="diagnosisuremia"
    )
    await check(
        result,
        "660fe42c7c13e47524344eaf5b7d2e7c",
        {"session_info": ["limma_3.58.1"], "results": ["1007_s_at\t-0.27752", "1053_at\t-0.0715477"]},
    )


@env.task
async def test_limma_array_exclude_samples() -> None:
    # upstream case: "test_limma_differential - exclude_samples"
    args = "--sample_id_col name --exclude_samples_col name --exclude_samples_values Dima_1;Dima_9"
    result = await run_array(args, **DIAGNOSIS)
    await check(
        result,
        ARRAY_MODEL,
        {"session_info": ["limma_3.58.1"], "results": ["1007_s_at\t-0.4778282", "1053_at\t-0.08958278"]},
    )


@env.task
async def test_limma_array_subset_to_contrast() -> None:
    # upstream case: "test_limma_differential - subset_to_contrast"
    result = await run_array("--sample_id_col name --subset_to_contrast_samples TRUE", **DIAGNOSIS)
    await check(
        result,
        ARRAY_MODEL,
        {"session_info": ["limma_3.58.1"], "results": ["1007_s_at\t-0.2775254", "1053_at\t-0.071547786"]},
    )


async def run_counts(args: str) -> LimmaDifferentialResult:
    return await limma_differential(
        await fixture(EXPRESSION + "SRP254919.samplesheet.csv"),
        await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv"),
        prefix="test",
        contrast_variable="treatment",
        reference="hND6",
        target="mCherry",
        args="--sample_id_col sample --probe_id_col gene_id --use_voom TRUE" + args,
    )


NORMALISED = ["ENSMUSG00000023978\t6.11247620232167", "ENSMUSG00000059991\t4.52751370160052"]


@env.task
async def test_limma_voom_blocking() -> None:
    # upstream case: "test_limma_differential - voom_blocking"
    result = await run_counts(" --blocking_variables sample_number")
    await check(
        result,
        "cb9d919594a83c8d15e3693bd599a0bb",
        {
            "results": ["ENSMUSG00000023978\t-4.89014922224", "ENSMUSG00000059991\t0.77922"],
            "normalised_counts": NORMALISED,
        },
    )


@env.task
async def test_limma_voom() -> None:
    # upstream case: "test_limma_differential - voom"
    result = await run_counts("")
    await check(
        result,
        SHEET_MODEL,
        {
            "results": ["ENSMUSG00000023978\t-2.84055986312942", "ENSMUSG00000059991\t0.7055"],
            "normalised_counts": NORMALISED,
        },
    )


@env.task
async def test_limma_voom_mixed() -> None:
    # upstream case: "test_limma_differential - voom_mixed"
    result = await run_counts(" --block sample_number")
    await check(
        result,
        SHEET_MODEL,
        {
            "results": ["ENSMUSG00000023978\t-2.8363444336503", "ENSMUSG00000059991\t0.7360"],
            "normalised_counts": NORMALISED,
        },
    )


@env.task
async def test_limma_voom_mixed_seed() -> None:
    # upstream case: "test_limma_differential - voom_mixed with seed"
    result = await run_counts(" --block sample_number --seed 1 --round_digits 2")
    await check(result, SHEET_MODEL)
    await assert_md5(result.results, "3cafd117ece830f7d5ef78c8f5813c48", label="results")
    assert result.normalised_counts is not None, "no normalised counts"
    await assert_md5(result.normalised_counts, "1c953ffa43ca6d882db45c7d8be86fd5", label="normalised counts")


tests = [
    test_limma_formula_interaction,
    test_limma_array,
    test_limma_array_formula,
    test_limma_array_exclude_samples,
    test_limma_array_subset_to_contrast,
    test_limma_voom_blocking,
    test_limma_voom,
    test_limma_voom_mixed,
    test_limma_voom_mixed_seed,
]
