"""Tests for flyte_bio.modules.blast.

The cases and expected md5s are upstream's blast/makeblastdb and blast/blastn
tests (modules commit 9a1cc9f, BLAST 2.17.0+). As upstream, the database's
``.nin`` and ``.njs`` (which record its creation date) are only checked for
presence, and blastn's report (which records the database's date) for its
query line.
"""

import hashlib

from flyte.io import Dir, File

from flyte_bio.modules.blast import blast_blastn, blast_makeblastdb
from tests.framework import env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"
DB_FILES = [f"genome.fasta.{ext}" for ext in ("ndb", "nhr", "nin", "njs", "not", "nsq", "ntf", "nto")]
STABLE_DB_MD5S = {
    "genome.fasta.ndb": "0d553c830656469211de113c5022f06d",
    "genome.fasta.nhr": "f4b4ddb034fd3dd7b25c89e9d50c004e",
    "genome.fasta.not": "1e53e9d08f1d23af0299cfa87478a7bb",
    "genome.fasta.nsq": "982cbc7d9e38743b9b1037588862b9da",
    "genome.fasta.ntf": "de1250813f0c7affc6d12dac9d0fb6bb",
    "genome.fasta.nto": "33cdeccccebe80329f1fdbee7f5874cb",
}
QUERY_LINE = "Query= MT192765.1 Severe acute respiratory syndrome coronavirus 2 isolate"


async def check_db(db: Dir, label: str) -> None:
    md5s = {}
    async for file in db.walk():
        async with file.open("rb") as fh:
            md5s[file.name] = hashlib.md5(bytes(await fh.read())).hexdigest()
    assert sorted(md5s) == DB_FILES, f"{label}: files {sorted(md5s)}, expected {DB_FILES}"
    for name, md5 in STABLE_DB_MD5S.items():
        assert md5s[name] == md5, f"{label}: {name} md5 {md5s[name]}, expected {md5}"


async def check_report(report: File, label: str) -> None:
    assert report.name == "test.txt", f"{label}: output named {report.name}"
    async with report.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    assert QUERY_LINE in text, f"{label}: report lacks {QUERY_LINE!r}"


@env.task
async def test_blast_makeblastdb_fasta() -> None:
    await check_db(await blast_makeblastdb(await fixture(GENOME), args="-dbtype nucl"), "fasta")


@env.task
async def test_blast_makeblastdb_gzipped_fasta() -> None:
    await check_db(await blast_makeblastdb(await fixture(GENOME + ".gz"), args="-dbtype nucl"), "fasta.gz")


@env.task
async def test_blast_blastn_fasta() -> None:
    db = await blast_makeblastdb(await fixture(GENOME), args="-dbtype nucl")
    await check_report(await blast_blastn(await fixture(GENOME), db, prefix="test"), "fasta")


@env.task
async def test_blast_blastn_gzipped_fasta() -> None:
    db = await blast_makeblastdb(await fixture(GENOME), args="-dbtype nucl")
    await check_report(await blast_blastn(await fixture(GENOME + ".gz"), db, prefix="test"), "fasta.gz")


tests = [
    test_blast_makeblastdb_fasta,
    test_blast_makeblastdb_gzipped_fasta,
    test_blast_blastn_fasta,
    test_blast_blastn_gzipped_fasta,
]
