"""Tests for flyte_bio.modules.rundbcan.

The cases and expected values are upstream's rundbcan/database and
rundbcan/cazymeannotation tests (modules commit 9a1cc9f, run_dbcan 5.2.9):

- the full database download, checking the database folder's entries, as
  upstream does;
- the annotation of a Bacteroides fragilis region's proteins against the
  test-datasets' mini database (``dbcan.tar.gz``), checking upstream's md5s,
  the overview's line count and a GH127 pectin row, and a GH76 HMM hit.
"""

import tarfile

import flyte
from flyte.io import File

from flyte_bio.modules.rundbcan import rundbcan_cazymeannotation, rundbcan_database
from tests.framework import assert_md5, env, fixture

DATABASE_ENTRIES = [
    "CAZy.dmnd",
    "PUL.dmnd",
    "STP.hmm",
    "TCDB.dmnd",
    "TF.dmnd",
    "TF.hmm",
    "dbCAN-PUL",
    "dbCAN-PUL.xlsx",
    "dbCAN-sub.hmm",
    "dbCAN.hmm",
    "fam-substrate-mapping.tsv",
    "peptidase_db.dmnd",
    "sulfatlas_db.dmnd",
]


@env.task
async def rundbcan_archive_entries(archive: File) -> list[str]:
    """The sorted entries of ``dbcan_db/`` in ``archive``, streamed."""
    entries = set()
    local = await archive.download()
    with tarfile.open(local, "r|gz") as tar:
        for member in tar:
            parts = member.name.removeprefix("./").split("/")
            if len(parts) >= 2 and parts[0] == "dbcan_db" and parts[1]:
                entries.add(parts[1])
    return sorted(entries)


@env.task
async def test_rundbcan_database() -> None:
    archive = await rundbcan_database()
    assert archive.name == "dbcan_db.tar.gz", f"database archive is {archive.name}"
    # The database is a few GiB: list it in a pod with room for it (downloading it
    # needs more than 2Gi of memory).
    lister = rundbcan_archive_entries.override(resources=flyte.Resources(cpu=1, memory="8Gi", disk="20Gi"))
    entries = await lister(archive=archive)
    assert entries == DATABASE_ENTRIES, f"dbcan_db holds {entries}"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_rundbcan_cazymeannotation_faa() -> None:
    result = await rundbcan_cazymeannotation(
        await fixture("genomics/prokaryotes/bacteroides_fragilis/genome/region_2500000-3000000.faa.gz"),
        await fixture("genomics/prokaryotes/bacteroides_fragilis/genome/db/dbcan.tar.gz"),
        prefix="test",
    )
    overview = (await text(result.cazyme_annotation)).splitlines()
    assert len(overview) == 58, f"overview has {len(overview)} lines, expected 58"
    assert any(line.startswith("WP_005784055.1\t") and line.endswith("\tGH127_e5\tpectin") for line in overview), (
        "overview lacks WP_005784055.1's GH127_e5 pectin row"
    )
    assert "GH76" in await text(result.dbcanhmm_results), "dbCAN HMM results lack GH76"
    await assert_md5(result.cazyme_annotation, "980562e2e1bdf7b21bba4d8b967ed434", label="test_overview.tsv")
    await assert_md5(result.dbcandiamond_results, "58c89af41c26aab794f6f5bcc44af6f1", label="test_diamond.out")
    await assert_md5(result.dbcanhmm_results, "e6ca4efe7353732ffa5ddb537c895d9a", label="test_dbCAN_hmm_results.tsv")
    await assert_md5(result.dbcansub_results, "2a5e25f3e9d25a8753b9de9103711448", label="test_dbCANsub_hmm_results.tsv")


tests = [
    test_rundbcan_database,
    test_rundbcan_cazymeannotation_faa,
]
