"""Tests for flyte_bio.modules.deseq2.

The cases and expected md5s are upstream's deseq2/differential tests (modules
commit 9a1cc9f), with each case's test config ``args``. The contrast-sheet
cases run both of upstream's contrasts (``treatment_mCherry_hND6_`` without
blocking, ``treatment_mCherry_hND6_sample_number`` blocking on
``sample_number``), as upstream's channel does; the "csv" case, like
upstream's single-item matrix channel, runs only the first. As upstream, the
plots and RDS are checked by name, and the session info for the DESeq2
version.
"""

import tempfile
from dataclasses import dataclass
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.deseq2 import Deseq2DifferentialResult, deseq2_differential
from tests.framework import assert_md5, env, fixture

DREAM = "https://github.com/nf-core/test-datasets/raw/differentialabundance/modules_testdata/variancepartition_dream/"
EXPRESSION = "genomics/mus_musculus/rnaseq_expression/"
VERSIONS_MD5 = "62c34d4f3d11c2bdbe2e4f7e6086b35e"
UNBLOCKED = "treatment_mCherry_hND6_"
BLOCKED = "treatment_mCherry_hND6_sample_number"
# Upstream's model.txt md5s for the two sheet contrasts, the same in every sheet case.
UNBLOCKED_MODEL = "d2113d82b76046c319e6602da2ad74d6"
BLOCKED_MODEL = "fa05126a58cb67c107d45426b0bdea83"


@dataclass
class Expected:
    results: str
    size_factors: str
    normalised_counts: str
    model: str
    rlog: str | None = None
    vst: str | None = None


async def check(result: Deseq2DifferentialResult, prefix: str, expected: Expected) -> None:
    await assert_md5(result.results, expected.results, label=f"{prefix} results")
    await assert_md5(result.size_factors, expected.size_factors, label=f"{prefix} size factors")
    await assert_md5(result.normalised_counts, expected.normalised_counts, label=f"{prefix} normalised counts")
    await assert_md5(result.model, expected.model, label=f"{prefix} model")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")
    for name, file, md5 in (("rlog", result.rlog_counts, expected.rlog), ("vst", result.vst_counts, expected.vst)):
        if md5 is None:
            assert file is None, f"{prefix}: unexpected {name} counts"
        else:
            assert file is not None, f"{prefix}: no {name} counts"
            await assert_md5(file, md5, label=f"{prefix} {name}")
    assert result.dispersion_plot_png.name == f"{prefix}.deseq2.dispersion.png"
    assert result.dispersion_plot_pdf.name == f"{prefix}.deseq2.dispersion.pdf"
    assert result.rdata.name == f"{prefix}.dds.rld.rds"
    async with result.session_info.open("rb") as fh:
        session = bytes(await fh.read()).decode()
    assert "DESeq2_1.34.0" in session, f"{prefix}: session info lacks DESeq2_1.34.0"


async def run_formula(prefix: str, formula: str, comparison: str, args: str) -> Deseq2DifferentialResult:
    return await deseq2_differential(
        await fixture(DREAM + "metadata.tsv"),
        await fixture(DREAM + "counts.tsv"),
        prefix=prefix,
        formula=formula,
        comparison=comparison,
        args=args,
    )


INTERACTION_ARGS = "--round_digits 2 --sample_id_col sample --vs_method rlog --shrink_lfc FALSE --seed 1"
SHRINK_ARGS = "--round_digits 2 --sample_id_col sample --vs_method rlog --shrink_lfc TRUE --seed 1"
APEGLM_ARGS = (
    "--round_digits 2 --sample_id_col sample --vs_method rlog --shrink_lfc TRUE --shrink_lfc_type apeglm --seed 1"
)
INTERACTION = "genotype_WT_KO_treatment_Control_Treated"
GENOTYPE = "genotype_WT_KO"
INTERACTION_SHARED = {
    "size_factors": "0367d1ef76c5b9f01fc1f89647809b5d",
    "normalised_counts": "6d9d7c02be9169ba4a649f8fb0f0b7ea",
    "rlog": "0e9d0074e4acc651f50234fcce55a1cf",
}


@env.task
async def test_deseq2_formula_comparison_interaction() -> None:
    # upstream case: "RNAseq - Feature Counts - formula + comparison contrast string - interaction"
    result = await run_formula(INTERACTION, "~ genotype * treatment", "genotypeWT.treatmentTreated", INTERACTION_ARGS)
    expected = Expected(
        "5612fc7e92a7a9b497b9198471ed3602", model="82165b8b24da8f8cffafa457290feeb3", **INTERACTION_SHARED
    )
    await check(result, INTERACTION, expected)


@env.task
async def test_deseq2_formula_limma_contrast() -> None:
    # upstream case: "RNAseq - Feature Counts - formula + limma contrast string - interaction"
    result = await run_formula(GENOTYPE, "~ 0 + genotype", "genotypeWT - genotypeKO", INTERACTION_ARGS)
    expected = Expected(
        "59256d9c5b4f12e4422dcf80776da841", model="227c742304821b46168b1ebf10adfdce", **INTERACTION_SHARED
    )
    await check(result, GENOTYPE, expected)


@env.task
async def test_deseq2_formula_limma_contrast_shrunken() -> None:
    # upstream case: "RNAseq - Feature Counts - formula + limma contrast string - shrunken"
    result = await run_formula(GENOTYPE, "~ 0 + genotype", "genotypeWT - genotypeKO", SHRINK_ARGS)
    expected = Expected(
        "23511305c222de275f888fee844f055c", model="227c742304821b46168b1ebf10adfdce", **INTERACTION_SHARED
    )
    await check(result, GENOTYPE, expected)


@env.task
async def test_deseq2_apeglm_shrinkage() -> None:
    # upstream case: "RNAseq - Feature Counts - apeglm shrinkage"
    result = await run_formula(INTERACTION, "~ genotype * treatment", "genotypeWT.treatmentTreated", APEGLM_ARGS)
    expected = Expected(
        "2e6f15c26629c232c7a052999eee8a49", model="82165b8b24da8f8cffafa457290feeb3", **INTERACTION_SHARED
    )
    await check(result, INTERACTION, expected)


async def expect_apeglm_error(**contrast: str) -> None:
    try:
        await deseq2_differential(
            await fixture(DREAM + "metadata.tsv"),
            await fixture(DREAM + "counts.tsv"),
            prefix=GENOTYPE,
            args=APEGLM_ARGS,
            **contrast,
        )
    except Exception as err:
        assert "apeglm shrinkage requires a model coefficient" in str(err), f"unexpected error: {err}"
    else:
        raise AssertionError("apeglm shrinkage without a model coefficient succeeded")


@env.task
async def test_deseq2_standard_contrast_apeglm_error() -> None:
    # upstream case: "RNAseq - Feature Counts - standard contrast - apeglm error"
    await expect_apeglm_error(contrast_variable="genotype", reference="KO", target="WT")


@env.task
async def test_deseq2_limma_contrast_apeglm_error() -> None:
    # upstream case: "RNAseq - Feature Counts - formula + limma contrast string - apeglm error"
    await expect_apeglm_error(formula="~ 0 + genotype", comparison="genotypeWT - genotypeKO")


async def run_sheet(
    args: str,
    expected: dict[str, Expected],
    counts: File | None = None,
    control_genes: File | None = None,
    transcript_lengths: File | None = None,
) -> None:
    """Run upstream's sheet contrasts in ``expected``; ``args`` takes each one's blocking variable as ``{}``."""
    sheet = await fixture(EXPRESSION + "SRP254919.samplesheet.csv")
    counts = counts or await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv")
    for prefix, blocking in ((UNBLOCKED, ""), (BLOCKED, "sample_number")):
        if prefix not in expected:
            continue
        result = await deseq2_differential(
            sheet,
            counts,
            prefix=prefix,
            contrast_variable="treatment",
            reference="mCherry",
            target="hND6",
            control_genes=control_genes,
            transcript_lengths=transcript_lengths,
            args=args.format(blocking),
        )
        await check(result, prefix, expected[prefix])


MATRIX_ARGS = "--round_digits 5 --blocking_variables {} --vs_method rlog"


def both(results: tuple[str, str], size_factors: str, normalised_counts: str, **counts: str) -> dict[str, Expected]:
    """Expected outputs of the two sheet contrasts, which share all but their results and model."""
    return {
        UNBLOCKED: Expected(results[0], size_factors, normalised_counts, UNBLOCKED_MODEL, **counts),
        BLOCKED: Expected(results[1], size_factors, normalised_counts, BLOCKED_MODEL, **counts),
    }


MATRIX_EXPECTED = both(
    ("791cdba2615a445cded13cae95df73ef", "2438053a4bdc869f467a12d3c22c7ba7"),
    "6332f21889ecac387bb12eeb04f23849",
    "46ab7200c626649ab6256ed797ef5071",
    rlog="b1adc1fba6bd0c8b55973608f4b97030",
)


@env.task
async def test_deseq2_contrasts_matrix() -> None:
    # upstream case: "mouse - contrasts - matrix"
    await run_sheet(MATRIX_ARGS, MATRIX_EXPECTED)


@env.task
async def test_deseq2_contrasts_matrix_lengths() -> None:
    # upstream case: "mouse - contrasts - matrix - lengths"
    lengths = await fixture(EXPRESSION + "SRP254919.spoofed_lengths.tsv")
    expected = both(
        ("944176b73455aa7c8de3ec32c03edef6", "6a9bc76c9d54034c90fa159372f97516"),
        "69183b94202a806067962d0df0b9875b",
        "7050f44c460cc13e3f101d048d503527",
        rlog="22a4b117246b2317e0f4daf7919703f2",
    )
    await run_sheet(MATRIX_ARGS, expected, transcript_lengths=lengths)


@env.task
async def test_deseq2_contrasts_matrix_no_blocking() -> None:
    # upstream case: "mouse - contrasts - matrix - no blocking". Upstream's config drops the blocking option
    # only for a null blocking value; the sheet's are empty strings or sample_number, so it is kept.
    await run_sheet("--blocking_variables {} --vs_method rlog --round_digits 5", MATRIX_EXPECTED)


@env.task
async def test_deseq2_contrasts_matrix_spikes() -> None:
    # upstream case: "mouse - contrasts - matrix - spikes"
    spikes = await fixture(EXPRESSION + "SRP254919.spikes.tsv")
    expected = both(
        ("51ba90cf073d9d5f38ebb043d7e5b220", "a04e2df3712d5be27e39c50f83319ec4"),
        "1708374146b4e8c28fbc78fa292102cd",
        "6f3af146a6f8fbb8f565dd66e79d2cb1",
        rlog="39df6fb2969afc009a602261e864238f",
    )
    args = "--round_digits 5 --sizefactors_from_controls TRUE --blocking_variables {} --vs_method rlog"
    await run_sheet(args, expected, control_genes=spikes)


@env.task
async def test_deseq2_contrasts_matrix_strip_spikes() -> None:
    # upstream case: "mouse - contrasts - matrix - strip spikes"
    spikes = await fixture(EXPRESSION + "SRP254919.spikes.tsv")
    expected = both(
        ("c4d269cace6dfe42e81e39d13b4b0354", "e1820951d5d5c40fa7fd6368e8f62e8c"),
        "a526a1bb5cb270f5a4a1641fefd289b3",
        "89c5922a30218a17ab6ee84e031f6a9b",
        rlog="4b4526ddcca010f5c75f841b2cef269f",
    )
    await run_sheet(MATRIX_ARGS, expected, control_genes=spikes)


@env.task
async def test_deseq2_contrasts_matrix_csv() -> None:
    # upstream case: "mouse - contrasts - matrix - csv". Upstream rewrites the matrix through splitCsv and
    # collectFile into test.tsv; its lines hold no commas, so that's the same content under a new name.
    original = await fixture(EXPRESSION + "SRP254919.salmon.merged.gene_counts.top1000cov.tsv")
    async with original.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert not any("," in line for line in lines), "the matrix holds commas: splitCsv would change it"
    path = Path(tempfile.mkdtemp()) / "test.tsv"
    path.write_text("".join(line + "\n" for line in lines))
    counts = await File.from_local(str(path))
    await run_sheet(MATRIX_ARGS, {UNBLOCKED: MATRIX_EXPECTED[UNBLOCKED]}, counts=counts)


@env.task
async def test_deseq2_contrasts_matrix_vst_nsub() -> None:
    # upstream case: "mouse - contrasts - matrix - vst nsub"
    expected = both(
        ("791cdba2615a445cded13cae95df73ef", "2438053a4bdc869f467a12d3c22c7ba7"),
        "6332f21889ecac387bb12eeb04f23849",
        "46ab7200c626649ab6256ed797ef5071",
        vst="41c39026998bb5943e19b6012e750d89",
    )
    await run_sheet("--round_digits 5 --blocking_variables {} --vs_method vst --vst_nsub 500", expected)


@env.task
async def test_deseq2_contrasts_matrix_subset_to_contrast() -> None:
    # upstream case: "mouse - contrasts - matrix - subset to contrast"
    args = "--round_digits 5 --blocking_variables {} --vs_method rlog --subset_to_contrast_samples TRUE"
    await run_sheet(args, MATRIX_EXPECTED)


@env.task
async def test_deseq2_contrasts_matrix_exclude_samples() -> None:
    # upstream case: "mouse - contrasts - matrix - exclude samples"
    expected = both(
        ("7d28d65b7ab94306c1328027a6531405", "f509c2de7d6935942eff594a0aca662e"),
        "87d6d634ffea2d7a9df9f9bae91b4d90",
        "90327e9bfc3bbafb96e11b13b3cd5e19",
        rlog="c500bcad14f845113f85610d1ebf643c",
    )
    args = (
        "--round_digits 4 --blocking_variables {} --vs_method rlog "
        "--exclude_samples_col sample_number --exclude_samples_values sample1"
    )
    await run_sheet(args, expected)


tests = [
    test_deseq2_formula_comparison_interaction,
    test_deseq2_formula_limma_contrast,
    test_deseq2_formula_limma_contrast_shrunken,
    test_deseq2_apeglm_shrinkage,
    test_deseq2_standard_contrast_apeglm_error,
    test_deseq2_limma_contrast_apeglm_error,
    test_deseq2_contrasts_matrix,
    test_deseq2_contrasts_matrix_lengths,
    test_deseq2_contrasts_matrix_no_blocking,
    test_deseq2_contrasts_matrix_spikes,
    test_deseq2_contrasts_matrix_strip_spikes,
    test_deseq2_contrasts_matrix_csv,
    test_deseq2_contrasts_matrix_vst_nsub,
    test_deseq2_contrasts_matrix_subset_to_contrast,
    test_deseq2_contrasts_matrix_exclude_samples,
]
