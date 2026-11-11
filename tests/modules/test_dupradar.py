"""Tests for flyte_bio.modules.dupradar.

The duplication matrix is checked against upstream's dupradar snapshot md5s
(the PDFs aren't byte-stable, so — as upstream — only their presence is
checked).

The two MultiQC tables each end in values that are zero to within rounding: an
intercept of about 2e-11, and a fitted duplication rate of about 2e-9. Their
trailing digits vary between node CPUs, so these tables are checked line by
line, with those values compared within a tolerance. The expected values are
upstream's: tables printing them hash to upstream's snapshot md5s.
"""

import hashlib
from dataclasses import dataclass

from flyte.io import File

from flyte_bio.modules.dupradar import DupradarResult, dupradar
from tests.framework import assert_md5, env, fixture

SARSCOV2 = "genomics/sarscov2/"
# As upstream's test config: the sarscov2 GTF has CDS but no exon features.
FEATURE_TYPE = "--feature_type CDS"

# md5 of the curve table's lines before its data line (identical in both cases).
CURVE_HEADER_MD5 = "dd8bcbc47917f2251f6b7edadde1dd60"


@dataclass
class Expected:
    matrix_md5: str
    intercept: float  # upstream's intercept table hashes to intercept_md5 with this value
    intercept_md5: str
    curve_x: float  # upstream's curve table hashes to curve_md5 with these values
    curve_y: float
    curve_md5: str


SINGLE_END = Expected(
    matrix_md5="2beda4c140548a2b8c91bf6bde01ddc6",
    intercept=2.14334458552868e-11,
    intercept_md5="602a30c09b79533abbbb76efbb168e3e",
    curve_x=1.04739460591778,
    curve_y=2.14334458548274e-09,
    curve_md5="b41f82f7d515d5e2b3e7993add90d26b",
)
PAIRED_END = Expected(
    matrix_md5="5a327feaba56e5ea96a2cb7e8577e196",
    intercept=2.14334458552858e-11,
    intercept_md5="4071a87d8db9731e9a74e51f00523056",
    curve_x=1.04739460591778,
    curve_y=2.14334458548265e-09,
    curve_md5="af58f18dd239447f04fa4a313124a549",
)


async def read_lines(table: File) -> list[str]:
    async with table.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def check_intercept_near(table: File, expected: float, label: str) -> None:
    """Check the MultiQC intercept table: header lines exact, value within 1e-9 of ``expected``."""
    lines = await read_lines(table)
    assert lines[0] == "#id: DupInt", f"{label}: unexpected header {lines[0]!r}"
    assert lines[-2] == "Sample dupRadar_intercept", f"{label}: unexpected column line {lines[-2]!r}"
    sample, value = lines[-1].split()
    assert sample == "test", f"{label}: unexpected sample {sample!r}"
    assert abs(float(value) - expected) < 1e-9, f"{label}: intercept {value}, expected about {expected}"


async def check_curve_near(table: File, x: float, y: float, label: str) -> None:
    """Check the MultiQC curve table: header lines by md5, the data point within tolerance."""
    lines = await read_lines(table)
    header = "\n".join(lines[:-1]) + "\n"
    got = hashlib.md5(header.encode()).hexdigest()
    assert got == CURVE_HEADER_MD5, f"{label}: header md5 {got}, expected {CURVE_HEADER_MD5}"
    gx, gy = (float(v) for v in lines[-1].split())
    assert abs(gx - x) <= 1e-9 * abs(x), f"{label}: expression {gx}, expected about {x}"
    assert abs(gy - y) < 1e-12, f"{label}: duplication rate {gy}, expected about {y}"


async def check(result: DupradarResult, expected: Expected, label: str) -> None:
    await assert_md5(result.dup_matrix, expected.matrix_md5, label=f"{label} dupMatrix")
    await check_intercept_near(result.intercept_mqc, expected.intercept, label=f"{label} dup_intercept_mqc")
    await check_curve_near(
        result.curve_mqc, expected.curve_x, expected.curve_y, label=f"{label} duprateExpDensCurve_mqc"
    )
    for name in ("test_duprateExpBoxplot.pdf", "test_expressionHist.pdf", "test_duprateExpDens.pdf"):
        assert await result.results.get_file(name) is not None, f"{label}: no {name}"


@env.task
async def test_dupradar_single_end() -> None:
    # upstream case: dupradar "sarscov2 - bam - single_end"
    bam = await fixture(SARSCOV2 + "illumina/bam/test.single_end.bam")
    gtf = await fixture(SARSCOV2 + "genome/genome.gtf")
    result = await dupradar(bam, gtf, prefix="test", strandedness="forward", single_end=True, args=FEATURE_TYPE)
    await check(result, SINGLE_END, "single_end")


@env.task
async def test_dupradar_paired_end() -> None:
    # upstream case: dupradar "sarscov2 - bam - paired_end"
    bam = await fixture(SARSCOV2 + "illumina/bam/test.paired_end.bam")
    gtf = await fixture(SARSCOV2 + "genome/genome.gtf")
    result = await dupradar(bam, gtf, prefix="test", strandedness="forward", single_end=False, args=FEATURE_TYPE)
    await check(result, PAIRED_END, "paired_end")


tests = [test_dupradar_single_end, test_dupradar_paired_end]
