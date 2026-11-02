"""Tests for flyte_bio.modules.gunzip."""

import asyncio
import hashlib

from flyte_bio.modules.gunzip import gunzip
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/"


@env.task
async def test_gunzip() -> None:
    # Decompressed, and named like the input without .gz (as upstream's gunzip).
    archive, plain = await asyncio.gather(fixture(GENOME + "genome.fasta.gz"), fixture(GENOME + "genome.fasta"))
    out = await gunzip(archive=archive)
    assert out.path.endswith("/genome.fasta"), out.path
    async with plain.open("rb") as fh:
        expected = hashlib.md5(bytes(await fh.read())).hexdigest()
    await assert_md5(out, expected, label="genome.fasta")


tests = [test_gunzip]
