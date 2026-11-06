"""Tests for flyte_bio.modules.picard.

Upstream's picard/markduplicates snapshot only records the output name and
the metrics header (which embeds run-specific paths), so this checks the
same shape: the marked BAM is produced and the metrics file is a Picard
duplication-metrics report for this command.

Upstream's picard/collectmultiplemetrics snapshots (modules commit 9a1cc9f,
Picard 3.5.0) record only the sorted names of the metrics files and PDFs,
since every metrics file's header records the command line and date. These
tests check the same names, plus that each metrics file is a Picard report of
this command.
"""

from flyte.io import File

from flyte_bio.modules.picard import CollectMultipleMetricsResult, picard_collectmultiplemetrics, picard_markduplicates
from tests.framework import assert_nonempty, env, fixture

SC2 = "genomics/sarscov2/"
COLLECTMULTIPLEMETRICS_METRICS = [
    "test.CollectMultipleMetrics.alignment_summary_metrics",
    "test.CollectMultipleMetrics.base_distribution_by_cycle_metrics",
    "test.CollectMultipleMetrics.insert_size_metrics",
    "test.CollectMultipleMetrics.quality_by_cycle_metrics",
    "test.CollectMultipleMetrics.quality_distribution_metrics",
]
COLLECTMULTIPLEMETRICS_PDFS = [
    "test.CollectMultipleMetrics.base_distribution_by_cycle.pdf",
    "test.CollectMultipleMetrics.insert_size_histogram.pdf",
    "test.CollectMultipleMetrics.quality_by_cycle.pdf",
    "test.CollectMultipleMetrics.quality_distribution.pdf",
    "test.CollectMultipleMetrics.read_length_histogram.pdf",
]

COMMAND = "# CollectMultipleMetrics --INPUT test.paired_end.sorted.bam "


async def text_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


@env.task
async def test_markduplicates_sorted_bam() -> None:
    # upstream case: picard/markduplicates "sarscov2 [sorted bam]"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam")
    result = await picard_markduplicates(bam, prefix="test.md", args="--ASSUME_SORT_ORDER queryname")
    await assert_nonempty(result.bam, label="marked bam")
    lines = await text_lines(result.metrics)
    assert lines[0] == "## htsjdk.samtools.metrics.StringHeader", lines[:2]
    assert lines[1].startswith("# MarkDuplicates --INPUT ") and "--OUTPUT " in lines[1], lines[1]
    assert any(line.startswith("## METRICS CLASS\tpicard.sam.DuplicationMetrics") for line in lines), "no metrics table"


async def check_collectmultiplemetrics(result: CollectMultipleMetricsResult, reference: bool) -> None:
    assert [f.name for f in result.metrics] == COLLECTMULTIPLEMETRICS_METRICS, [f.name for f in result.metrics]
    assert [f.name for f in result.pdfs] == COLLECTMULTIPLEMETRICS_PDFS, [f.name for f in result.pdfs]
    for metrics in result.metrics:
        lines = await text_lines(metrics)
        assert lines[0] == "## htsjdk.samtools.metrics.StringHeader", (metrics.name, lines[:2])
        command = lines[1]
        assert command.startswith(COMMAND), (metrics.name, command)
        assert ("--REFERENCE_SEQUENCE genome.fasta" in command) == reference, (metrics.name, command)
    for pdf in result.pdfs:
        async with pdf.open("rb") as fh:
            assert bytes(await fh.read(5)) == b"%PDF-", f"{pdf.name} is not a PDF"


async def sorted_bam() -> tuple[File, File]:
    bam = await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam")
    return bam, await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai")


@env.task
async def test_picard_collectmultiplemetrics() -> None:
    # upstream case: picard/collectmultiplemetrics "test-picard-collectmultiplemetrics"
    bam, bai = await sorted_bam()
    fasta = await fixture(SC2 + "genome/genome.fasta")
    result = await picard_collectmultiplemetrics(bam, prefix="test", bai=bai, fasta=fasta)
    await check_collectmultiplemetrics(result, reference=True)


@env.task
async def test_picard_collectmultiplemetrics_nofasta() -> None:
    # upstream case: picard/collectmultiplemetrics "test-picard-collectmultiplemetrics-nofasta"
    bam, bai = await sorted_bam()
    result = await picard_collectmultiplemetrics(bam, prefix="test", bai=bai)
    await check_collectmultiplemetrics(result, reference=False)


@env.task
async def test_picard_collectmultiplemetrics_cram() -> None:
    # upstream case: picard/collectmultiplemetrics "test-picard-collectmultiplemetrics-cram". Despite its
    # name, upstream passes the BAM, with the sarscov2 FASTA and a homo_sapiens .fai beside it.
    bam, bai = await sorted_bam()
    fasta = await fixture(SC2 + "genome/genome.fasta")
    fai = await fixture("genomics/homo_sapiens/genome/genome.fasta.fai")
    result = await picard_collectmultiplemetrics(bam, prefix="test", bai=bai, fasta=fasta, fai=fai)
    await check_collectmultiplemetrics(result, reference=True)


tests = [
    test_markduplicates_sorted_bam,
    test_picard_collectmultiplemetrics,
    test_picard_collectmultiplemetrics_nofasta,
    test_picard_collectmultiplemetrics_cram,
]
