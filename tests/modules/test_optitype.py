"""Tests for flyte_bio.modules.optitype.

The cases and expected values are upstream's optitype tests (modules commit
9a1cc9f, OptiType 1.5.0): the HLA example BAM typed with GLPK and with CBC,
both with ``--enumerate 1 --beta 0.009``. As upstream, the snapshot is the
calls' B and C alleles and read count (upstream skips locus A, whose calls
aren't stable on this data) and the coverage plot's name.
"""

import asyncio

from flyte.io import File

from flyte_bio.modules.optitype import OptitypeResult, optitype
from tests.framework import env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/example_hla_pe.sorted.bam"
ARGS = "--enumerate 1 --beta 0.009"
EXPECTED_CALLS = ["B*08:01", "B*57:01", "C*06:02", "C*07:01", "1439.0"]


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check_calls(result: OptitypeResult, label: str) -> None:
    assert result.hla_type.name == "test_result.tsv", f"{label}: calls are {result.hla_type.name}"
    assert result.coverage_plot.name == "test_coverage_plot.pdf", f"{label}: plot is {result.coverage_plot.name}"
    # upstream: readLines()[1].split('\t')[3..7]
    calls = (await text(result.hla_type)).splitlines()[1].split("\t")[3:8]
    assert calls == EXPECTED_CALLS, f"{label}: calls {calls}, expected {EXPECTED_CALLS}"


async def run_optitype(solver: str) -> OptitypeResult:
    bam, bai = await asyncio.gather(fixture(BAM), fixture(BAM + ".bai"))
    return await optitype([bam], [bai], "test", seq_type="dna", solver=solver, args=ARGS)


@env.task
async def test_optitype_glpk() -> None:
    # upstream case: optitype "hla - bam"
    await check_calls(await run_optitype("glpk"), "glpk")


@env.task
async def test_optitype_cbc() -> None:
    # upstream case: optitype "hla - bam - cbc"
    await check_calls(await run_optitype("cbc"), "cbc")


tests = [test_optitype_glpk, test_optitype_cbc]
