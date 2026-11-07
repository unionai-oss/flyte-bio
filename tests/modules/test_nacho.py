"""Tests for flyte_bio.modules.nacho.

The cases and expected md5s are upstream's nacho/qc and nacho/normalize tests
(modules commit 9a1cc9f): the four salmon RCC files of the nanostring test
data, with their sample sheet. As upstream's snapshots: the QC's report and
plot names (their bytes aren't reproducible) and its two tables' md5s; the
GEO-normalised counts' md5s.
"""

import csv

from flyte.io import File

from flyte_bio.modules.nacho import nacho_normalize, nacho_qc
from tests.framework import assert_md5, env, fixture

SAMPLE_SHEET = "https://raw.githubusercontent.com/nf-core/test-datasets/nanostring/samplesheets/samplesheet_test.csv"
PLOTS = [
    "AVG_vs_BD_mqc.png",
    "AVG_vs_MED_mqc.png",
    "BD_mqc.png",
    "FOV_mqc.png",
    "HKF_mqc.png",
    "HK_mqc.png",
    "LOD_mqc.png",
    "Neg_mqc.png",
    "PCA1_vs_PCA2_mqc.png",
    "PCA_mqc.png",
    "PCAi_mqc.png",
    "POSF_vs_NEGF_mqc.png",
    "Pos_mqc.png",
    "Pos_vs_neg_mqc.png",
    "Posctrl_linearity_mqc.png",
    "plot_normf_mqc.png",
]


async def rcc_files() -> tuple[list[File], File]:
    """The RCC files the sample sheet's RCC_FILE column points at, and the sheet itself."""
    sheet = await fixture(SAMPLE_SHEET)
    async with sheet.open("rb") as fh:
        rows = list(csv.DictReader(bytes(await fh.read()).decode().splitlines()))
    return [await fixture(row["RCC_FILE"]) for row in rows], sheet


@env.task
async def test_nacho_qc_salmon() -> None:
    rccs, sheet = await rcc_files()
    result = await nacho_qc(rccs, sheet)
    assert [f.name for f in result.reports] == ["NanoQC.html", "NanoQC_with_outliers.html"], result.reports
    assert sorted(f.name for f in result.plots) == sorted(PLOTS), [f.name for f in result.plots]
    tables = {f.name: f for f in result.tables}
    assert sorted(tables) == ["hk_detected_mqc.txt", "normalized_qc_mqc.txt"], sorted(tables)
    await assert_md5(tables["hk_detected_mqc.txt"], "61209383acc2abf6fc3ea309b5a5e094", label="hk_detected_mqc.txt")
    await assert_md5(tables["normalized_qc_mqc.txt"], "e9f9c9342d7abf167b02a0147cfcd9be", label="normalized_qc_mqc.txt")


@env.task
async def test_nacho_normalize_salmon() -> None:
    rccs, sheet = await rcc_files()
    result = await nacho_normalize(rccs, sheet, prefix="test", args="--norm_method GEO")
    await assert_md5(result.normalized_counts, "4e0aa3cab4d7206e432dc7fbb0c7b334", label="test.tsv")
    await assert_md5(result.normalized_counts_wo_hk, "4ca3ffb0e0b15ca8b9dbe553714cb199", label="test_wo_HKnorm.tsv")


tests = [test_nacho_qc_salmon, test_nacho_normalize_salmon]
