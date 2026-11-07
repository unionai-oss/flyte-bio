"""Tests for flyte_bio.modules.nanocomp.

The cases are upstream's nanocomp tests (modules commit 9a1cc9f, NanoComp
1.25.6). As upstream, the HTML reports and ``NanoStats.txt`` aren't
reproducible, so the tests check the outputs' names against upstream's
snapshot, and that ``NanoStats.txt`` has more than one line. NanoComp also
writes two weighted histograms that upstream declares no output for (so its
snapshot doesn't list them); the module keeps them, and the tests expect them. The input-type
checks upstream makes in its script are tested here too.
"""

from flyte.io import File

from flyte_bio.modules.nanocomp import NanocompResult, file_name, input_type, nanocomp
from tests.framework import env, fixture

COMMON_HTML = [
    "test_NanoComp-report.html",
    "test_NanoComp_N50.html",
    "test_NanoComp_OverlayHistogram.html",
    "test_NanoComp_OverlayHistogram_Normalized.html",
    "test_NanoComp_OverlayHistogram_Weighted.html",  # not an upstream output
    "test_NanoComp_OverlayLogHistogram.html",
    "test_NanoComp_OverlayLogHistogram_Normalized.html",
    "test_NanoComp_OverlayLogHistogram_Weighted.html",  # not an upstream output
    "test_NanoComp_lengths_violin.html",
    "test_NanoComp_log_length_violin.html",
    "test_NanoComp_number_of_reads.html",
    "test_NanoComp_quals_violin.html",
    "test_NanoComp_total_throughput.html",
]


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: NanocompResult, html: list[str]) -> None:
    names = [file_name(f) for f in result.html]
    assert names == sorted(html), f"NanoComp html: {names}, expected {sorted(html)}"
    assert file_name(result.report_html) == "test_NanoComp-report.html", file_name(result.report_html)
    assert file_name(result.stats_txt) == "test_NanoStats.txt", file_name(result.stats_txt)
    lines = len((await text(result.stats_txt)).splitlines())
    assert lines > 1, f"test_NanoStats.txt has {lines} lines"


@env.task
async def test_nanocomp_fastq() -> None:
    reads = await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz")
    await check(await nanocomp([reads], prefix="test"), COMMON_HTML)


@env.task
async def test_nanocomp_summary() -> None:
    summary = await fixture("genomics/sarscov2/nanopore/sequencing_summary/test.sequencing_summary.txt")
    html = [
        *COMMON_HTML,
        "test_NanoComp_ActivePoresOverTime.html",
        "test_NanoComp_CumulativeYieldPlot_Gigabases.html",
        "test_NanoComp_sequencing_speed_over_time.html",
    ]
    await check(await nanocomp([summary], prefix="test"), html)


@env.task
async def test_nanocomp_input_types() -> None:
    assert input_type("test.fastq.gz") == "fastq"
    assert input_type("test.fa") == "fasta"
    assert input_type("test.sorted.bam") == "bam"
    assert input_type("test.sequencing_summary.txt") == "summary"
    for name in ("test", "test.vcf.gz"):
        try:
            input_type(name)
        except ValueError:
            pass
        else:
            raise AssertionError(f"nanocomp accepted {name}")


tests = [test_nanocomp_fastq, test_nanocomp_summary, test_nanocomp_input_types]
