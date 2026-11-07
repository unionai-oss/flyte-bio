"""Tests for flyte_bio.modules.toulligqc.

The cases and expectations are upstream's toulligqc tests (modules commit
9a1cc9f, ToulligQC 2.9.1): the report is named ``report.html``, and lines
12-75 of ``report.data`` (12-68 for the FASTQ case) match upstream's
snapshot. The earlier lines record the run's date, paths and versions. As
upstream's test config, the barcoded case passes ``--barcoding`` and
barcodes 01-12.
"""

import os

from flyte.io import File

from flyte_bio.modules.toulligqc import toulligqc
from tests.framework import env, fixture

NANOPORE = "genomics/sarscov2/nanopore/"
SUMMARY = NANOPORE + "sequencing_summary/test2.sequencing_summary.txt"
BARCODES = ",".join(f"barcode{i:02d}" for i in range(1, 13))

# Upstream's snapshot lines, with this common key prefix factored out.
B = "basecaller.sequencing.summary.1d.extractor."

SUMMARY_LINES = [
    "sequencing.telemetry.extractor.software.analysis=1d_basecalling",
    B + "read.count=100",
    B + "read.pass.count=100",
    B + "read.fail.count=0",
    B + "read.pass.ratio=1.0",
    B + "read.fail.ratio=0.0",
    B + "read.count.frequency=100",
    B + "read.pass.frequency=100.0",
    B + "read.fail.frequency=0.0",
    B + "yield=38253",
    B + "n50=365",
    B + "l50=67",
    B + "run.time=29006.4915",
    B + "channel.occupancy.statistics.count=85.0",
    B + "channel.occupancy.statistics.mean=1.1764705882352942",
    B + "channel.occupancy.statistics.std=0.4412231534591759",
    B + "channel.occupancy.statistics.min=1.0",
    B + "channel.occupancy.statistics.25%=1.0",
    B + "channel.occupancy.statistics.50%=1.0",
    B + "channel.occupancy.statistics.75%=1.0",
    B + "channel.occupancy.statistics.max=3.0",
    B + "all.read.length.count=100.0",
    B + "all.read.length.mean=382.53",
    B + "all.read.length.std=219.3747982818722",
    B + "all.read.length.min=228.0",
    B + "all.read.length.25%=271.75",
    B + "all.read.length.50%=305.5",
    B + "all.read.length.75%=418.0",
    B + "all.read.length.max=1664.0",
    B + "pass.reads.sequence.length.mean=382.53",
    B + "pass.reads.sequence.length.std=219.3747982818722",
    B + "pass.reads.sequence.length.min=228.0",
    B + "pass.reads.sequence.length.25%=271.75",
    B + "pass.reads.sequence.length.50%=305.5",
    B + "pass.reads.sequence.length.75%=418.0",
    B + "pass.reads.sequence.length.max=1664.0",
    B + "fail.reads.sequence.length.mean=nan",
    B + "fail.reads.sequence.length.std=nan",
    B + "fail.reads.sequence.length.min=nan",
    B + "fail.reads.sequence.length.25%=nan",
    B + "fail.reads.sequence.length.50%=nan",
    B + "fail.reads.sequence.length.75%=nan",
    B + "fail.reads.sequence.length.max=nan",
    B + "all.read.qscore.mean=12.084663391113281",
    B + "all.read.qscore.std=1.7714887857437134",
    B + "all.read.qscore.min=7.7340922355651855",
    B + "all.read.qscore.25%=11.089608192443848",
    B + "all.read.qscore.50%=11.957954406738281",
    B + "all.read.qscore.75%=13.238139390945435",
    B + "all.read.qscore.max=17.272123336791992",
    B + "pass.reads.mean.qscore.mean=12.084663391113281",
    B + "pass.reads.mean.qscore.std=1.7714887857437134",
    B + "pass.reads.mean.qscore.min=7.7340922355651855",
    B + "pass.reads.mean.qscore.25%=11.089608192443848",
    B + "pass.reads.mean.qscore.50%=11.957954406738281",
    B + "pass.reads.mean.qscore.75%=13.238139390945435",
    B + "pass.reads.mean.qscore.max=17.272123336791992",
    B + "fail.reads.mean.qscore.mean=nan",
    B + "fail.reads.mean.qscore.std=nan",
    B + "fail.reads.mean.qscore.min=nan",
    B + "fail.reads.mean.qscore.25%=nan",
    B + "fail.reads.mean.qscore.50%=nan",
    B + "fail.reads.mean.qscore.75%=nan",
    B + "fail.reads.mean.qscore.max=nan",
]

SUMMARY_BARCODES_LINES = [
    "toulligqc.info.extractor.duration=0.0",
    "sequencing.telemetry.extractor.software.analysis=1d_basecalling",
    B + "read.count=100",
    B + "read.pass.count=100",
    B + "read.fail.count=0",
    B + "read.pass.ratio=1.0",
    B + "read.fail.ratio=0.0",
    B + "read.count.frequency=100",
    B + "read.pass.frequency=100.0",
    B + "read.fail.frequency=0.0",
    B + "yield=38253",
    B + "n50=365",
    B + "l50=67",
    B + "run.time=29006.4915",
    B + "channel.occupancy.statistics.count=85.0",
    B + "channel.occupancy.statistics.mean=1.1764705882352942",
    B + "channel.occupancy.statistics.std=0.4412231534591759",
    B + "channel.occupancy.statistics.min=1.0",
    B + "channel.occupancy.statistics.25%=1.0",
    B + "channel.occupancy.statistics.50%=1.0",
    B + "channel.occupancy.statistics.75%=1.0",
    B + "channel.occupancy.statistics.max=3.0",
    B + "all.read.length.count=100.0",
    B + "all.read.length.mean=382.53",
    B + "all.read.length.std=219.3747982818722",
    B + "all.read.length.min=228.0",
    B + "all.read.length.25%=271.75",
    B + "all.read.length.50%=305.5",
    B + "all.read.length.75%=418.0",
    B + "all.read.length.max=1664.0",
    B + "pass.reads.sequence.length.mean=382.53",
    B + "pass.reads.sequence.length.std=219.3747982818722",
    B + "pass.reads.sequence.length.min=228.0",
    B + "pass.reads.sequence.length.25%=271.75",
    B + "pass.reads.sequence.length.50%=305.5",
    B + "pass.reads.sequence.length.75%=418.0",
    B + "pass.reads.sequence.length.max=1664.0",
    B + "fail.reads.sequence.length.mean=nan",
    B + "fail.reads.sequence.length.std=nan",
    B + "fail.reads.sequence.length.min=nan",
    B + "fail.reads.sequence.length.25%=nan",
    B + "fail.reads.sequence.length.50%=nan",
    B + "fail.reads.sequence.length.75%=nan",
    B + "fail.reads.sequence.length.max=nan",
    B + "all.read.qscore.mean=12.084663391113281",
    B + "all.read.qscore.std=1.7714887857437134",
    B + "all.read.qscore.min=7.7340922355651855",
    B + "all.read.qscore.25%=11.089608192443848",
    B + "all.read.qscore.50%=11.957954406738281",
    B + "all.read.qscore.75%=13.238139390945435",
    B + "all.read.qscore.max=17.272123336791992",
    B + "pass.reads.mean.qscore.mean=12.084663391113281",
    B + "pass.reads.mean.qscore.std=1.7714887857437134",
    B + "pass.reads.mean.qscore.min=7.7340922355651855",
    B + "pass.reads.mean.qscore.25%=11.089608192443848",
    B + "pass.reads.mean.qscore.50%=11.957954406738281",
    B + "pass.reads.mean.qscore.75%=13.238139390945435",
    B + "pass.reads.mean.qscore.max=17.272123336791992",
    B + "fail.reads.mean.qscore.mean=nan",
    B + "fail.reads.mean.qscore.std=nan",
    B + "fail.reads.mean.qscore.min=nan",
    B + "fail.reads.mean.qscore.25%=nan",
    B + "fail.reads.mean.qscore.50%=nan",
    B + "fail.reads.mean.qscore.75%=nan",
]

FASTQ_LINES = [
    "sequencing.telemetry.extractor.run.id=97ca62ca093ff43533aa34c38a10b1d6325e7e7b",
    "sequencing.telemetry.extractor.sample.id=RN20097",
    "sequencing.telemetry.extractor.model.file=Unknown",
    "sequencing.telemetry.extractor.flowcell.id=FAP51364",
    B + "read.count=100",
    B + "read.pass.count=96",
    B + "read.fail.count=4",
    B + "read.pass.ratio=0.96",
    B + "read.fail.ratio=0.04",
    B + "read.count.frequency=100",
    B + "read.pass.frequency=96.0",
    B + "read.fail.frequency=4.0",
    B + "yield=38253",
    B + "n50=365",
    B + "l50=67",
    B + "run.time=2035.0",
    B + "channel.occupancy.statistics.count=85.0",
    B + "channel.occupancy.statistics.mean=1.1764705882352942",
    B + "channel.occupancy.statistics.std=0.4412231534591759",
    B + "channel.occupancy.statistics.min=1.0",
    B + "channel.occupancy.statistics.25%=1.0",
    B + "channel.occupancy.statistics.50%=1.0",
    B + "channel.occupancy.statistics.75%=1.0",
    B + "channel.occupancy.statistics.max=3.0",
    B + "all.read.length.count=100.0",
    B + "all.read.length.mean=382.53",
    B + "all.read.length.std=219.37479828187222",
    B + "all.read.length.min=228.0",
    B + "all.read.length.25%=271.75",
    B + "all.read.length.50%=305.5",
    B + "all.read.length.75%=418.0",
    B + "all.read.length.max=1664.0",
    B + "pass.reads.sequence.length.mean=385.6145833333333",
    B + "pass.reads.sequence.length.std=222.9568931568171",
    B + "pass.reads.sequence.length.min=228.0",
    B + "pass.reads.sequence.length.25%=272.0",
    B + "pass.reads.sequence.length.50%=312.0",
    B + "pass.reads.sequence.length.75%=418.0",
    B + "pass.reads.sequence.length.max=1664.0",
    B + "fail.reads.sequence.length.mean=308.5",
    B + "fail.reads.sequence.length.std=79.87698875980074",
    B + "fail.reads.sequence.length.min=261.0",
    B + "fail.reads.sequence.length.25%=267.75",
    B + "fail.reads.sequence.length.50%=272.5",
    B + "fail.reads.sequence.length.75%=313.25",
    B + "fail.reads.sequence.length.max=428.0",
    B + "all.read.qscore.mean=12.82710075378418",
    B + "all.read.qscore.std=1.9102991819381714",
    B + "all.read.qscore.min=8.199999809265137",
    B + "all.read.qscore.25%=11.702500104904175",
    B + "all.read.qscore.50%=12.714999675750732",
    B + "all.read.qscore.75%=14.082499980926514",
    B + "all.read.qscore.max=18.110000610351562",
    B + "pass.reads.mean.qscore.mean=13.002917289733887",
    B + "pass.reads.mean.qscore.std=1.737709403038025",
    B + "pass.reads.mean.qscore.min=9.84000015258789",
    B + "pass.reads.mean.qscore.25%=11.90749979019165",
]

BAM_LINES = [
    "sequencing.telemetry.extractor.run.id=Unknown",
    "sequencing.telemetry.extractor.sample.id=Unknown",
    "sequencing.telemetry.extractor.model.file=Unknown",
    "sequencing.telemetry.extractor.software.name=minimap2",
    "sequencing.telemetry.extractor.software.version=2.17-r974-dirty",
    "sequencing.telemetry.extractor.flowcell.id=Unknown",
    "sequencing.telemetry.extractor.basecalling.date=Unknown",
    "sequencing.telemetry.extractor.pass.threshold.qscore=9",
    B + "read.count=100",
    B + "read.pass.count=96",
    B + "read.fail.count=4",
    B + "read.pass.ratio=0.96",
    B + "read.fail.ratio=0.04",
    B + "read.count.frequency=100",
    B + "read.pass.frequency=96.0",
    B + "read.fail.frequency=4.0",
    B + "yield=38253",
    B + "n50=365",
    B + "l50=67",
    B + "run.time=99.0",
    B + "channel.occupancy.statistics.count=1.0",
    B + "channel.occupancy.statistics.mean=100.0",
    B + "channel.occupancy.statistics.std=nan",
    B + "channel.occupancy.statistics.min=100.0",
    B + "channel.occupancy.statistics.25%=100.0",
    B + "channel.occupancy.statistics.50%=100.0",
    B + "channel.occupancy.statistics.75%=100.0",
    B + "channel.occupancy.statistics.max=100.0",
    B + "all.read.length.count=100.0",
    B + "all.read.length.mean=382.53",
    B + "all.read.length.std=219.37479828187222",
    B + "all.read.length.min=228.0",
    B + "all.read.length.25%=271.75",
    B + "all.read.length.50%=305.5",
    B + "all.read.length.75%=418.0",
    B + "all.read.length.max=1664.0",
    B + "pass.reads.sequence.length.mean=385.6145833333333",
    B + "pass.reads.sequence.length.std=222.9568931568171",
    B + "pass.reads.sequence.length.min=228.0",
    B + "pass.reads.sequence.length.25%=272.0",
    B + "pass.reads.sequence.length.50%=312.0",
    B + "pass.reads.sequence.length.75%=418.0",
    B + "pass.reads.sequence.length.max=1664.0",
    B + "fail.reads.sequence.length.mean=308.5",
    B + "fail.reads.sequence.length.std=79.87698875980074",
    B + "fail.reads.sequence.length.min=261.0",
    B + "fail.reads.sequence.length.25%=267.75",
    B + "fail.reads.sequence.length.50%=272.5",
    B + "fail.reads.sequence.length.75%=313.25",
    B + "fail.reads.sequence.length.max=428.0",
    B + "all.read.qscore.mean=12.82710075378418",
    B + "all.read.qscore.std=1.9102991819381714",
    B + "all.read.qscore.min=8.199999809265137",
    B + "all.read.qscore.25%=11.702500104904175",
    B + "all.read.qscore.50%=12.714999675750732",
    B + "all.read.qscore.75%=14.082499980926514",
    B + "all.read.qscore.max=18.110000610351562",
    B + "pass.reads.mean.qscore.mean=13.002917289733887",
    B + "pass.reads.mean.qscore.std=1.737709403038025",
    B + "pass.reads.mean.qscore.min=9.84000015258789",
    B + "pass.reads.mean.qscore.25%=11.90749979019165",
    B + "pass.reads.mean.qscore.50%=12.820000171661377",
    B + "pass.reads.mean.qscore.75%=14.09500002861023",
    B + "pass.reads.mean.qscore.max=18.110000610351562",
]


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def check(ont_file: File, expected: list[str], args: str = "") -> None:
    result = await toulligqc(ont_file, prefix="test", args=args)
    assert result.report_html is not None, "toulligqc wrote no report html"
    assert file_name(result.report_html) == "report.html", f"report html named {file_name(result.report_html)}"
    async with result.report_data.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    actual = lines[11 : 11 + len(expected)]
    for i, (got, want) in enumerate(zip(actual, expected, strict=False)):
        assert got == want, f"report.data line {12 + i}: expected {want!r}, got {got!r}"
    assert len(actual) == len(expected), f"report.data has {len(lines)} lines, expected at least {11 + len(expected)}"


@env.task
async def test_toulligqc_sequencing_summary() -> None:
    await check(await fixture(SUMMARY), SUMMARY_LINES)


@env.task
async def test_toulligqc_sequencing_summary_barcodes() -> None:
    await check(await fixture(SUMMARY), SUMMARY_BARCODES_LINES, args=f"--barcoding --barcodes {BARCODES}")


@env.task
async def test_toulligqc_fastq() -> None:
    await check(await fixture(NANOPORE + "fastq/test.fastq.gz"), FASTQ_LINES)


@env.task
async def test_toulligqc_bam() -> None:
    await check(await fixture(NANOPORE + "bam/test.sorted.bam"), BAM_LINES)


tests = [
    test_toulligqc_sequencing_summary,
    test_toulligqc_sequencing_summary_barcodes,
    test_toulligqc_fastq,
    test_toulligqc_bam,
]
