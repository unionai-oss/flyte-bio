"""Tests for flyte_bio.modules.picard.

Upstream's picard/markduplicates snapshot only records the output name and
the metrics header (which embeds run-specific paths), so this checks the
same shape: the marked BAM is produced and the metrics file is a Picard
duplication-metrics report for this command.
"""

from flyte_bio.modules.picard import picard_markduplicates
from tests.framework import assert_nonempty, env, fixture


@env.task
async def test_markduplicates_sorted_bam() -> None:
    # upstream case: picard/markduplicates "sarscov2 [sorted bam]"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam")
    result = await picard_markduplicates(bam, prefix="test.md", args="--ASSUME_SORT_ORDER queryname")
    await assert_nonempty(result.bam, label="marked bam")
    async with result.metrics.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert lines[0] == "## htsjdk.samtools.metrics.StringHeader", lines[:2]
    assert lines[1].startswith("# MarkDuplicates --INPUT ") and "--OUTPUT " in lines[1], lines[1]
    assert any(line.startswith("## METRICS CLASS\tpicard.sam.DuplicationMetrics") for line in lines), "no metrics table"


tests = [test_markduplicates_sorted_bam]
