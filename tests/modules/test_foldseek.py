"""Tests for flyte_bio.modules.foldseek.

The cases are upstream's foldseek/createdb and foldseek/easysearch tests
(modules commit 9a1cc9f, Foldseek 10.941cd33). createdb's expected md5s are
upstream's snapshot. easysearch searches 8tim against a database of 1tim built
by createdb (named ``test_db``, as upstream's setup); upstream treats its
``.m8`` as unstable and checks it holds the expected hit, as these tests do.
The second case names a database that isn't there, so the search falls back
to the single database in the Dir.
"""

import hashlib

from flyte.io import Dir

from flyte_bio.modules.foldseek import foldseek_createdb, foldseek_easysearch
from tests.framework import env, fixture

QUERY = "proteomics/pdb/8tim.pdb"
TARGET = "proteomics/pdb/1tim.pdb"
EXPECTED_HIT = "8tim_A\t1tim_A\t0.967\t247\t8\t0\t1\t247\t1\t247\t1.152E-43\t1526"
CREATEDB = {
    "test": "5509f0851612be656f0659d62e290366",
    "test.dbtype": "f1d3ff8443297732862df21dc4e57262",
    "test.index": "e7eb45c3cbd57553ff8f024c7c83ed11",
    "test.lookup": "3e653ea1ae44ba5404342b37c72f05f6",
    "test.source": "1144e80754b5ba7bd85670e2be1c96db",
    "test_ca": "a5c5ffa9e813f7c35ed9a2737d4ef786",
    "test_ca.dbtype": "3fd85f9ee7ca8882c8caa747d0eef0b3",
    "test_ca.index": "4020bfeefff484c9c973a895fa76c4bd",
    "test_h": "eb1f0b213589e93fa639d3ae91bcff2f",
    "test_h.dbtype": "740bab4f9ec8808aedb68d6b1281aeb2",
    "test_h.index": "11940a92d3235f31b8a487c59a1ff812",
    "test_ss": "88d06580a2ee6a6a81a7821a423013a9",
    "test_ss.dbtype": "f1d3ff8443297732862df21dc4e57262",
    "test_ss.index": "e7eb45c3cbd57553ff8f024c7c83ed11",
}


async def target_database() -> Dir:
    return await foldseek_createdb(await fixture(TARGET), prefix="test_db")


async def check_search(db_name: str) -> None:
    aln = await foldseek_easysearch(
        await fixture(QUERY), await target_database(), db_name=db_name, prefix="test_search"
    )
    assert aln.name == "test_search.m8", f"output is {aln.name}"
    async with aln.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert EXPECTED_HIT in lines, f"test_search.m8 lacks {EXPECTED_HIT!r}: {lines[:5]}"


@env.task
async def test_foldseek_createdb_pdb() -> None:
    db = await foldseek_createdb(await fixture(TARGET), prefix="test")
    actual = {}
    async for file in db.walk():
        async with file.open("rb") as fh:
            actual[file.name] = hashlib.md5(bytes(await fh.read())).hexdigest()
    assert sorted(actual) == sorted(CREATEDB), f"files {sorted(actual)}"
    for name, md5 in CREATEDB.items():
        assert actual[name] == md5, f"{name} md5 {actual[name]}, expected {md5}"


@env.task
async def test_foldseek_easysearch_pdb() -> None:
    await check_search("test_db")


@env.task
async def test_foldseek_easysearch_lookup_fallback() -> None:
    await check_search("other_db")


tests = [
    test_foldseek_createdb_pdb,
    test_foldseek_easysearch_pdb,
    test_foldseek_easysearch_lookup_fallback,
]
