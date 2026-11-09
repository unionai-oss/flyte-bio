"""Tests for flyte_bio.modules.checkm2.

The cases and expected md5s are upstream's checkm2 tests (modules commit
9a1cc9f, CheckM2 1.1.0). As upstream:

- prediction runs with ``--lowmem`` on the Bacteroides fragilis genome
  against the reduced test database (a plain ``.dmnd`` File); the module's
  4 CPUs match upstream's test CPUs;
- the database download fetches the full version-3 database from Zenodo
  (about 1.7 GB) and checks only that it was written.
"""

from flyte_bio.modules.checkm2 import checkm2_databasedownload, checkm2_predict
from tests.framework import assert_md5, assert_nonempty, env, fixture


@env.task
async def test_checkm2_databasedownload() -> None:
    db = await checkm2_databasedownload()
    assert db.version == "3", f"database version {db.version}, expected 3"
    assert db.database.name == "checkm2_db_v3.dmnd", f"database named {db.database.name}"
    await assert_nonempty(db.database, label="database")


@env.task
async def test_checkm2_predict() -> None:
    result = await checkm2_predict(
        await fixture("genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz"),
        db=await fixture("delete_me/binette_checkm2/binette_checkm2_small_db.dmnd"),
        prefix="test",
        args="--lowmem",
    )
    await assert_md5(result.report, "d126e04913601c553bd1128735f31345", label="checkm2 report")


tests = [test_checkm2_databasedownload, test_checkm2_predict]
