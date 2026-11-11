"""Tests for flyte_bio.modules.decoupler.

The cases and checks are upstream's decoupler/decoupler tests (modules commit
9a1cc9f), with the test config's prefix (``deseq2_<id>``): the output file
names, the activity tables' headers and leading values, and the tool
versions upstream's snapshot records.
"""

from flyte.io import File

from flyte_bio.modules.decoupler import DecouplerResult, decoupler
from tests.framework import env, fixture

DIFFERENTIAL = "https://github.com/nf-core/test-datasets/raw/differentialabundance/modules_testdata/"
VERSIONS = (
    '"DECOUPLER_DECOUPLER":\n'
    "    decoupler-py: 1.8.0\n"
    "    matplotlib: 3.10.1\n"
    "    pandas: 2.2.3\n"
    "    python: 3.12.9\n"
    "    scanpy: 1.11.0\n"
)
PROGENY_HEADER = "\tAndrogen\tEGFR\tEstrogen\tHypoxia\tJAK-STAT\tMAPK\tNFkB\tPI3K\tTGFb\tTNFa\tTrail\tVEGF\tWNT\tp53"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: DecouplerResult, prefix: str, method: str, header: str, expected: dict[str, str]) -> None:
    names = [f"{prefix}_consensus_{{}}_decoupler.tsv", f"{prefix}_{method}_{{}}_decoupler.tsv"]
    assert [f.name for f in result.estimates] == [n.format("estimate") for n in names], result.estimates
    assert [f.name for f in result.pvals] == [n.format("pvals") for n in names], result.pvals
    plots = [f"{prefix}_consensus_estimate_decoupler_plot.png", f"{prefix}_{method}_estimate_decoupler_plot.png"]
    assert [f.name for f in result.plots] == plots, result.plots
    for table in [*result.estimates, *result.pvals]:
        assert header in await text(table), f"{table.name} lacks the header {header!r}"
    for table in [*result.estimates, *result.pvals]:
        for suffix, value in expected.items():
            if table.name.endswith(suffix):
                assert value in await text(table), f"{table.name} lacks {value!r}"
    assert await text(result.versions) == VERSIONS, await text(result.versions)


@env.task
async def test_decoupler_mouse_differential() -> None:
    # upstream case: "differential abundance data - normalized counts - mouse"
    result = await decoupler(
        await fixture(DIFFERENTIAL + "treatment_mCherry_hND6_.deseq2.results.tsv"),
        await fixture(DIFFERENTIAL + "progeny/mouse_network.tsv"),
        prefix="deseq2_treatment_mCherry_hND6",
        sample="treatment_mCherry_hND6",
        annot=await fixture(DIFFERENTIAL + "Mus_musculus.anno.feature_metadata.tsv"),
        args="--min_n 2 --transpose TRUE --ensembl_ids TRUE --methods mlm "
        "--features_id_col gene_id --features_symbol_col gene_name",
    )
    await check(
        result,
        "deseq2_treatment_mCherry_hND6",
        "mlm",
        PROGENY_HEADER,
        {
            "mlm_estimate_decoupler.tsv": "treatment_mCherry_hND6\t-0.0",
            "consensus_estimate_decoupler.tsv": "treatment_mCherry_hND6\t-0.0",
            "mlm_pvals_decoupler.tsv": "treatment_mCherry_hND6\t0.9",
            "consensus_pvals_decoupler.tsv": "treatment_mCherry_hND6\t0.9",
        },
    )


@env.task
async def test_decoupler_generic() -> None:
    # upstream case: "decoupler_test"
    result = await decoupler(
        await fixture("generic/tsv/expression.tsv"),
        await fixture("generic/tsv/network.tsv"),
        prefix="deseq2_test",
        sample="test",
        annot=await fixture(DIFFERENTIAL + "Mus_musculus.anno.feature_metadata.tsv"),
        args="--min_n 1 --features_id_col gene_id --features_symbol_col gene_name",
    )
    await check(
        result,
        "deseq2_test",
        "ulm",
        "\tT1\tT2\tT3\tT4\tT5",
        {
            "ulm_estimate_decoupler.tsv": "S01\t4.0",
            "consensus_estimate_decoupler.tsv": "S01\t1.1",
            "ulm_pvals_decoupler.tsv": "S01\t0.0",
            "consensus_pvals_decoupler.tsv": "S01\t0.2",
        },
    )


tests = [test_decoupler_mouse_differential, test_decoupler_generic]
