"""Tests for flyte_bio.modules.pangolin.

The cases are upstream's pangolin/updatedata and pangolin/run tests (modules
commit 9a1cc9f, pangolin 4.4 with pangolin-data 1.39 in the image). As
upstream: the run on the image's data checks the report's md5; updatedata
(live data) checks only that it ran, and the run on updated data checks the
report's name and header.
"""

from flyte_bio.modules.pangolin import pangolin_run, pangolin_updatedata
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"
HEADER = (
    "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,"
    "pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note"
)


@env.task
async def test_pangolin_run() -> None:
    # upstream case: pangolin/run "sarscov2 genome [fasta]"
    report = await pangolin_run(await fixture(GENOME), prefix="test")
    assert report.name == "test.pangolin.csv", report.name
    await assert_md5(report, "bdcdb7bcec59c2e6a12c11ac7ce43b09", label="test.pangolin.csv")


@env.task
async def test_pangolin_updatedata() -> None:
    # upstream case: pangolin/updatedata "str - testdb"
    await pangolin_updatedata()


@env.task
async def test_pangolin_run_updated_data() -> None:
    # upstream case: pangolin/run "sarscov2 genome [fasta] - update data"
    report = await pangolin_run(await fixture(GENOME), db=await pangolin_updatedata(), prefix="test")
    assert report.name == "test.pangolin.csv", report.name
    async with report.open("rb") as fh:
        first = bytes(await fh.read()).decode().splitlines()[0]
    assert HEADER in first, f"report header: {first}"


tests = [test_pangolin_run, test_pangolin_updatedata, test_pangolin_run_updated_data]
