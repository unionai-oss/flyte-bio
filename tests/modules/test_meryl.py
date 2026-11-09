"""Tests for flyte_bio.modules.meryl.

The cases and expected md5s are upstream's meryl/count, meryl/unionsum and
meryl/histogram tests (modules commit 9a1cc9f, meryl 1.4.1), with their
count steps chained in as upstream sets them up.

Upstream snapshots a database as its 129 files' ``name:md5,<md5>`` entries.
The tests here compare the md5 of those entries, sorted and joined with
newlines, computed from upstream's snapshot (the entries are listed in name
order there).
"""

import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.meryl import meryl_count, meryl_histogram, meryl_unionsum
from tests.framework import assert_md5, env, fixture

SARSCOV2 = "genomics/sarscov2/illumina/fastq/"
FRAGILIS = "genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_1.fastq.gz"


async def meryl_db_digest(db: File, name: str) -> str:
    """The md5 of the database's sorted ``<file>:md5,<md5>`` entries; ``db`` must hold ``name``."""
    assert db.name == f"{name}.tar.gz", f"database is {db.name}, expected {name}.tar.gz"
    root = Path(tempfile.mkdtemp())
    with tarfile.open(await db.download()) as tar:
        tar.extractall(root, filter="data")
    folder = root / name
    assert folder.is_dir(), f"{db.name} holds no {name}/"
    entries = sorted(f"{f.name}:md5,{hashlib.md5(f.read_bytes()).hexdigest()}" for f in folder.iterdir())
    assert len(entries) == 129, f"{name} has {len(entries)} files, expected 129"
    return hashlib.md5("\n".join(entries).encode()).hexdigest()


async def assert_meryl_db(db: File, name: str, expected: str) -> None:
    actual = await meryl_db_digest(db, name)
    assert actual == expected, f"{name}: database md5 mismatch (expected {expected}, got {actual})"


@env.task
async def test_meryl_count_fastq() -> None:
    dbs = await meryl_count([await fixture(FRAGILIS)], kvalue=21, prefix="test")
    assert len(dbs) == 1, f"meryl count wrote {len(dbs)} databases"
    await assert_meryl_db(dbs[0], "test.test1_1.meryl", "e4e0304b9349f176e69cfb3b187c1a02")


@env.task
async def test_meryl_unionsum_single_end() -> None:
    dbs = await meryl_count([await fixture(SARSCOV2 + "test_1.fastq.gz")], kvalue=21, prefix="test")
    db = await meryl_unionsum(dbs, kvalue=21, prefix="test")
    await assert_meryl_db(db, "test.unionsum.meryl", "1b21bd0289d99b1bdec54f6bab58d3e0")


@env.task
async def test_meryl_unionsum_paired_end() -> None:
    reads = [await fixture(SARSCOV2 + "test_1.fastq.gz"), await fixture(SARSCOV2 + "test_2.fastq.gz")]
    dbs = await meryl_count(reads, kvalue=21, prefix="test")
    assert [d.name for d in dbs] == ["test.test_1.meryl.tar.gz", "test.test_2.meryl.tar.gz"], [d.name for d in dbs]
    db = await meryl_unionsum(dbs, kvalue=21, prefix="test")
    await assert_meryl_db(db, "test.unionsum.meryl", "b0a63c8a6b2c1378973937540d5d1c95")


@env.task
async def test_meryl_histogram_fastq() -> None:
    dbs = await meryl_count([await fixture(FRAGILIS)], kvalue=21, prefix="test")
    hist = await meryl_histogram(dbs[0], kvalue=21, prefix="test")
    assert hist.name == "test.hist", f"meryl histogram wrote {hist.name}"
    await assert_md5(hist, "f75362ab9cd70d96621b3690e952085f", label="test.hist")


tests = [
    test_meryl_count_fastq,
    test_meryl_unionsum_single_end,
    test_meryl_unionsum_paired_end,
    test_meryl_histogram_fastq,
]
