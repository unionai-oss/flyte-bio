"""Tests for flyte_bio.modules.zip.

The case is upstream's zip test (modules commit 9a1cc9f). Upstream asserts
only that the archive exists (zip headers record times); this also checks
that it holds the input, byte for byte.
"""

import hashlib
import io
import zipfile

from flyte_bio.modules.zip import zip_files
from tests.framework import env, fixture

FASTA = "genomics/sarscov2/genome/genome.fasta"


@env.task
async def test_zip_fasta() -> None:
    # upstream case: "sarscov2 - fasta"
    fasta = await fixture(FASTA)
    archive = await zip_files([fasta], prefix="test")
    assert archive.name == "test.zip", f"unexpected archive name {archive.name}"
    async with archive.open("rb") as fh:
        data = bytes(await fh.read())
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = zf.namelist()
        members = [n for n in names if n.endswith("genome.fasta")]
        assert members == ["genome.fasta"], f"unexpected members {names}"
        zipped = zf.read(members[0])
    async with fasta.open("rb") as fh:
        original = bytes(await fh.read())
    assert hashlib.md5(zipped).hexdigest() == hashlib.md5(original).hexdigest(), "zipped FASTA differs"


tests = [test_zip_fasta]
