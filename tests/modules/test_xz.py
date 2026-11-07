"""Tests for flyte_bio.modules.xz.

The case and expected md5 are upstream's xz/decompress test (modules commit
9a1cc9f, xz 5.2.6). Upstream's setup step compresses the FASTA with its
xz/compress module (``xz -T 1 --stdout``); here the test compresses it with
Python's ``lzma`` (the same .xz container format, preset 6), since only the
decompressed bytes are checked.
"""

import asyncio
import lzma
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.xz import xz_decompress
from tests.framework import assert_md5, env, fixture


async def xz_compressed(file: File) -> File:
    """``file`` xz-compressed, as ``<name>.xz``."""
    async with file.open("rb") as fh:
        data = bytes(await fh.read())
    out = Path(tempfile.mkdtemp()) / f"{file.name or Path(file.path).name}.xz"
    await asyncio.to_thread(out.write_bytes, lzma.compress(data, format=lzma.FORMAT_XZ))
    return await File.from_local(str(out))


@env.task
async def test_xz_decompress_fasta() -> None:
    # upstream case: "sarscov2 - genome - fasta"
    archive = await xz_compressed(await fixture("genomics/sarscov2/genome/genome.fasta"))
    out = await xz_decompress(archive)
    assert out.path.endswith("/genome.fasta"), f"decompressed file is {out.path}, expected genome.fasta"
    await assert_md5(out, "6e9fe4042a72f2345f644f239272b7e6", label="genome.fasta")


tests = [test_xz_decompress_fasta]
