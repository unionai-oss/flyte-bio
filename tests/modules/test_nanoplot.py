"""Tests for flyte_bio.modules.nanoplot.

The cases are upstream's nanoplot tests (modules commit 9a1cc9f, NanoPlot
1.48.0). As upstream, the HTML reports aren't reproducible, so the tests
check their names against upstream's snapshot, and that ``NanoStats.txt``
has more than one line.
"""

from flyte.io import File

from flyte_bio.modules.nanoplot import NanoplotResult, file_name, nanoplot
from tests.framework import env, fixture


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: NanoplotResult, html: list[str]) -> None:
    names = [file_name(f) for f in result.html]
    assert names == html, f"NanoPlot html: {names}, expected {html}"
    assert [file_name(f) for f in result.txt] == ["NanoStats.txt"], f"NanoPlot txt: {result.txt}"
    lines = len((await text(result.txt[0])).splitlines())
    assert lines > 1, f"NanoStats.txt has {lines} lines"


@env.task
async def test_nanoplot_summary() -> None:
    summary = await fixture("genomics/sarscov2/nanopore/sequencing_summary/test.sequencing_summary.txt")
    html = [
        "ActivePores_Over_Time.html",
        "ActivityMap_ReadsPerChannel.html",
        "CumulativeYieldPlot_Gigabases.html",
        "CumulativeYieldPlot_NumberOfReads.html",
        "NanoPlot-report.html",
        "Non_weightedHistogramReadlength.html",
        "Non_weightedLogTransformed_HistogramReadlength.html",
        "NumberOfReads_Over_Time.html",
        "TimeLengthViolinPlot.html",
        "TimeQualityViolinPlot.html",
        "TimeSequencingSpeed_ViolinPlot.html",
        "WeightedHistogramReadlength.html",
        "WeightedLogTransformed_HistogramReadlength.html",
        "Yield_By_Length.html",
    ]
    await check(await nanoplot(summary), html)


@env.task
async def test_nanoplot_fastq() -> None:
    reads = await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz")
    html = [
        "LengthvsQualityScatterPlot_dot.html",
        "LengthvsQualityScatterPlot_kde.html",
        "NanoPlot-report.html",
        "Non_weightedHistogramReadlength.html",
        "Non_weightedLogTransformed_HistogramReadlength.html",
        "WeightedHistogramReadlength.html",
        "WeightedLogTransformed_HistogramReadlength.html",
        "Yield_By_Length.html",
    ]
    await check(await nanoplot(reads), html)


tests = [test_nanoplot_summary, test_nanoplot_fastq]
