"""Tests for flyte_bio.modules.mtmalign.

The cases are upstream's mtmalign/align tests (modules commit 9a1cc9f,
mTM-align 20220104) on the five seatoxin structures. mTM-align may be
nondeterministic, so, as upstream, they check that every structure is in the
alignment (and that the superposed structures were written).
"""

import gzip

from flyte.io import File

from flyte_bio.modules.mtmalign import mtmalign_align
from tests.framework import assert_nonempty, env, fixture_dir

STRUCTURES = (
    "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/structures/"
    "seatoxin-ref.tar.gz"
)
SEATOXIN_IDS = ("1ahl", "1apf", "1atx", "1bds", "1sh1")


async def seatoxin_pdbs() -> list[File]:
    structures = await fixture_dir(STRUCTURES, "seatoxin-ref")
    return sorted([f async for f in structures.walk()], key=lambda f: f.name)


async def check_alignment(file: File, name: str) -> None:
    assert file.name == name, f"output named {file.name}, expected {name}"
    async with file.open("rb") as fh:
        raw = bytes(await fh.read())
    text = (gzip.decompress(raw) if name.endswith(".gz") else raw).decode()
    missing = [i for i in SEATOXIN_IDS if f">{i}" not in text]
    assert not missing, f"{name} lacks {missing}"


@env.task
async def test_mtmalign_seatoxin_uncompressed() -> None:
    result = await mtmalign_align(await seatoxin_pdbs(), prefix="test")
    await check_alignment(result.alignment, "test.aln")
    assert result.structure.name == "test.pdb", f"structure named {result.structure.name}"
    await assert_nonempty(result.structure, label="test.pdb")


@env.task
async def test_mtmalign_seatoxin_compressed() -> None:
    result = await mtmalign_align(await seatoxin_pdbs(), prefix="test", compress=True)
    await check_alignment(result.alignment, "test.aln.gz")
    assert result.structure.name == "test.pdb.gz", f"structure named {result.structure.name}"


tests = [test_mtmalign_seatoxin_uncompressed, test_mtmalign_seatoxin_compressed]
