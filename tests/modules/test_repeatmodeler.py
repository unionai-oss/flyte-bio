"""Tests for flyte_bio.modules.repeatmodeler.

The cases and md5s are upstream's repeatmodeler/builddatabase and
repeatmodeler/repeatmodeler tests (modules commit 9a1cc9f, RepeatModeler
2.0.5). As upstream, the database's ``.nin`` and ``.njs`` (which record the
build date) are checked by presence, and RepeatModeler's alignments and log by
content.
"""

import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.repeatmodeler import repeatmodeler, repeatmodeler_builddatabase
from tests.framework import assert_md5, env, fixture

DATABASE_MD5S = {
    "test.nhr": "1a41cb6d0b00c28f62ad60e75ae2f6fc",
    "test.nin": None,
    "test.njs": None,
    "test.nnd": "2002e13acf59079a1a5782c918894579",
    "test.nni": "26a954ba0fd80983b550d8f6b8b35ff8",
    "test.nog": "30896f123998e926ea2237b89091e7fe",
    "test.nsq": "982cbc7d9e38743b9b1037588862b9da",
    "test.translation": "ccbb119522c09daa976a9015ba999329",
}


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


def member_md5s(archive: str) -> dict[str, str]:
    """Each file of a ``.tar.gz``, by name, with its md5."""
    root = Path(tempfile.mkdtemp())
    with tarfile.open(archive) as tar:
        tar.extractall(root, filter="data")
    return {p.name: hashlib.md5(p.read_bytes()).hexdigest() for p in root.iterdir()}


@env.task
async def test_repeatmodeler_builddatabase_sarscov2() -> None:
    # upstream case: repeatmodeler/builddatabase "sarscov2-genome_fasta"
    archive = await repeatmodeler_builddatabase(await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test")
    md5s = member_md5s(str(await archive.download()))
    assert sorted(md5s) == sorted(DATABASE_MD5S), sorted(md5s)
    for name, md5 in DATABASE_MD5S.items():
        assert md5 is None or md5s[name] == md5, f"{name}: md5 mismatch (expected {md5}, got {md5s[name]})"


@env.task
async def test_repeatmodeler_homo_sapiens() -> None:
    # upstream case: repeatmodeler/repeatmodeler "homo_sapiens-genome_fasta", from BuildDatabase as upstream
    archive = await repeatmodeler_builddatabase(
        await fixture("genomics/homo_sapiens/genome/genome.fasta"), prefix="test"
    )
    result = await repeatmodeler(archive, prefix="test")
    await assert_md5(result.fasta, "e25326771341204e1f8054d9529411e5", label="test.fa")
    assert "1 families discovered." in await text(result.log), "RepeatModeler log lacks '1 families discovered.'"
    assert (await text(result.stk)).startswith("# STOCKHOLM"), "alignments aren't Stockholm"


tests = [test_repeatmodeler_builddatabase_sarscov2, test_repeatmodeler_homo_sapiens]
