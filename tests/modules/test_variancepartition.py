"""Tests for flyte_bio.modules.variancepartition.

The cases and expected values are upstream's variancepartition/dream tests
(modules commit 9a1cc9f), with their ``args``: the model (and, where upstream
snapshots them, results and normalised counts) by md5, and the result
values upstream checks for. The contrast-sheet case runs the sheet's first
contrast only, as upstream's single-item matrix channel does. The module's
4 CPUs match the 4 that upstream's tests cap process_high at.
"""

from flyte.io import File

from flyte_bio.modules.variancepartition import DreamResult, variancepartition_dream
from tests.framework import assert_md5, env, fixture

DREAM = "https://github.com/nf-core/test-datasets/raw/differentialabundance/modules_testdata/variancepartition_dream/"
EXPRESSION = "genomics/mus_musculus/rnaseq_expression/"
VERSIONS_MD5 = "fc1f26eb2194018e99fc2916332676b7"
HEADER = "gene_id\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB"
INTERACTION = "genotype_WT_KO_treatment_Control_Treated"
INTERACTION_MODEL = "e67d5a9c42dd9374b370b6426832d08f"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: DreamResult, model: str, values: list[str]) -> None:
    await assert_md5(result.model, model, label="model")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")
    results = await text(result.results)
    for value in [HEADER, *values]:
        assert value in results, f"results lack {value!r}"


async def run_interaction(args: str) -> DreamResult:
    return await variancepartition_dream(
        await fixture(DREAM + "metadata.tsv"),
        await fixture(DREAM + "counts.tsv"),
        prefix=INTERACTION,
        formula="~ genotype * treatment",
        comparison="genotypeWT - treatmentTreated",
        args=args,
    )


async def run_sheet(prefix: str, **contrast: str) -> DreamResult:
    return await variancepartition_dream(
        await fixture(EXPRESSION + "SRP254919.samplesheet.csv"),
        await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv"),
        prefix=prefix,
        args="--round_digits 2",
        **contrast,
    )


@env.task
async def test_variancepartition_dream_interaction() -> None:
    # upstream case: "RNAseq - Feature Counts - formula + comparison contrast string - interaction"
    result = await run_interaction("--round_digits 2")
    await check(result, INTERACTION_MODEL, ["1309.83", "557\t279.88\t18.81"])


@env.task
async def test_variancepartition_dream_voom_interaction() -> None:
    # upstream case: "RNAseq - Voom - Feature Counts - formula + comparison contrast string - interaction"
    result = await run_interaction("--round_digits 2 --apply_voom true")
    await check(result, INTERACTION_MODEL, ["5.07", "2.97\t3.88\t38.17"])
    assert result.normalised_counts is not None, "no normalised counts"
    await assert_md5(result.normalised_counts, "84917b0fe7db41763ae3dc7c09e2ec4f", label="normalised counts")


@env.task
async def test_variancepartition_dream_formula_comparison() -> None:
    # upstream case: "Mus musculus - expression table - contrasts + formula + comparison contrast string"
    result = await run_sheet("treatment_mCherry_hND6", formula="~ treatment", comparison="treatmentmCherry")
    await check(result, "8cb944a29e3d3d0e540048b4d8331b17", ["-849.67", "-1050\t549\t-3.78"])


@env.task
async def test_variancepartition_dream_contrast_sheet() -> None:
    # upstream case: "Mus musculus - contrasts - matrix - no formula"
    result = await run_sheet(
        "treatment_mCherry_hND6_", contrast_variable="treatment", reference="mCherry", target="hND6"
    )
    await check(result, "30f74fcf5ad34ef1cc5505b0026963d6", ["1050", "849.67\t432.83\t3.8"])


@env.task
async def test_variancepartition_dream_weighted_comparison() -> None:
    # upstream case: "Mus musculus - expression table - contrasts + formula + weighted comparison contrast string"
    result = await run_sheet(
        "treatment_mCherry_hND6", formula="~ 0 + treatment", comparison="2 * treatmenthND6 - treatmentmCherry"
    )
    await check(result, "30f74fcf5ad34ef1cc5505b0026963d6", ["2124\t549\t4.84", "1707.33"])


@env.task
async def test_variancepartition_dream_no_intercept() -> None:
    # upstream case: "Mus musculus - expression table - contrasts + formula + comparison contrast string - no intercept"
    result = await run_sheet(
        "treatment_mCherry_hND6",
        formula="~ 0 + treatment + sample_number",
        comparison="treatmenthND6 - treatmentmCherry",
    )
    await check(result, "46fe9413749f1475edc6985647b0e4e1", ["-95.67", "1050\t549\t4.16"])


@env.task
async def test_variancepartition_dream_blocking_without_comparison() -> None:
    # upstream case: "Mus musculus - expression table - contrasts + blocking factors" (expected to fail)
    try:
        await run_sheet("treatment_mCherry_hND6", formula="~ treatment + (1 | sample_number)")
    except Exception as err:
        assert "Missing required options: contrast_string" in str(err), f"unexpected error: {err}"
    else:
        raise AssertionError("dream ran a formula without a comparison")


@env.task
async def test_variancepartition_dream_interaction_seed() -> None:
    # upstream case: "RNAseq - Feature Counts - formula + comparison contrast string - interaction with seed"
    result = await run_interaction("--seed 2 --round_digits 2")
    await assert_md5(result.results, "1d4c862325e4f5228c6277344c2f6b7d", label="results")
    await assert_md5(result.model, INTERACTION_MODEL, label="model")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")
    assert result.normalised_counts is None, "normalised counts without voom"


tests = [
    test_variancepartition_dream_interaction,
    test_variancepartition_dream_voom_interaction,
    test_variancepartition_dream_formula_comparison,
    test_variancepartition_dream_contrast_sheet,
    test_variancepartition_dream_weighted_comparison,
    test_variancepartition_dream_no_intercept,
    test_variancepartition_dream_blocking_without_comparison,
    test_variancepartition_dream_interaction_seed,
]
