"""Tests for flyte_bio.modules.mageck.

The cases and expected md5s are upstream's mageck count, mle and test tests
(modules commit 9a1cc9f, MAGeCK 0.5.9.5). Upstream leaves mle's gene summary
out of its snapshot (unstable), checking only that its header names genes.
"""

from flyte.io import File

from flyte_bio.modules.mageck import mageck_count, mageck_mle, mageck_test
from tests.framework import assert_md5, env, fixture

MAGECK = "genomics/mus_musculus/mageck/"


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def check_named(files: list[File], expected: dict[str, str], label: str) -> None:
    names = [file_name(f) for f in files]
    assert names == sorted(expected), f"{label}: wrote {names}, expected {sorted(expected)}"
    for file in files:
        await assert_md5(file, expected[file_name(file)], label=f"{label} {file_name(file)}")


@env.task
async def test_mageck_count_fastq() -> None:
    reads = [await fixture(MAGECK + "ERR376998.small.fastq.gz"), await fixture(MAGECK + "ERR376999.small.fastq.gz")]
    library = await fixture(MAGECK + "yusa_library.csv")
    result = await mageck_count(reads, library, sample_label="test,test2", prefix="test,test2")
    await assert_md5(result.count, "2759bdc010998f863518cf7a2df812dc", label="count")
    await assert_md5(result.norm, "d79f3863168ee0f0750a59d94f70ce1a", label="count_normalized")


@env.task
async def test_mageck_mle_csv() -> None:
    counts = await fixture(MAGECK + "count_table.csv")
    design = await fixture(MAGECK + "design_matrix.txt")
    result = await mageck_mle(counts, design, prefix="test")
    async with result.gene_summary[0].open("rb") as fh:
        header = bytes(await fh.read()).decode().splitlines()[0]
    assert "Gene" in header, f"gene summary header: {header!r}"
    await check_named(result.sgrna_summary, {"test.sgrna_summary.txt": "2d9c7b762ca4086f6d84a869f0961f3b"}, "sgrna")


@env.task
async def test_mageck_test_csv() -> None:
    result = await mageck_test(await fixture(MAGECK + "count_table.csv"), prefix="test", args="-t HL60.initial")
    await check_named(result.gene_summary, {"test.gene_summary.txt": "23d39164e44e7125a5f5912463fa5970"}, "gene")
    await check_named(result.sgrna_summary, {"test.sgrna_summary.txt": "26e4a50272a0289b0261a5380276f273"}, "sgrna")
    assert result.r_script is not None, "mageck test wrote no R script"
    await assert_md5(result.r_script, "9c913f55aa7e7142399396088332ea3b", label="R script")


@env.task
async def test_mageck_test_day0_csv() -> None:
    result = await mageck_test(
        await fixture(MAGECK + "count_table.csv"), prefix="test", args="--day0-label HL60.initial"
    )
    await check_named(
        result.gene_summary,
        {
            "test.HL60.final_vs_HL60.initial.gene_summary.txt": "973d57a2a432eb6c7e602ca466299b6d",
            "test.KBM7.final_vs_HL60.initial.gene_summary.txt": "5a39c0132560cf6a3818b0e463df99a9",
            "test.KBM7.initial_vs_HL60.initial.gene_summary.txt": "d532a4280945d28bad2186287b32cb3e",
            "test.gene_summary.txt": "2e5447f5930f01f9adea57bf89cec5fe",
        },
        "gene",
    )
    await check_named(
        result.sgrna_summary,
        {
            "test.HL60.final_vs_HL60.initial.sgrna_summary.txt": "a299a0f2cf7ce62bcf8aa4b4f48ea9db",
            "test.KBM7.final_vs_HL60.initial.sgrna_summary.txt": "aeb64ee78422ef4906161cb03c0c0ad2",
            "test.KBM7.initial_vs_HL60.initial.sgrna_summary.txt": "60c9104cad7d49bcead4d29cafb2af40",
            "test.sgrna_summary.txt": "60c9104cad7d49bcead4d29cafb2af40",
        },
        "sgrna",
    )
    assert result.r_script is not None, "mageck test wrote no R script"
    await assert_md5(result.r_script, "0849644a8ebfef2cf35d1ac351a02bca", label="R script")


tests = [test_mageck_count_fastq, test_mageck_mle_csv, test_mageck_test_csv, test_mageck_test_day0_csv]
