"""Tests for flyte_bio.modules.magus.

The cases are upstream's magus/align and magus/guidetree tests (modules commit
9a1cc9f, MAGUS 0.2.0) on the seatoxin sequences. MAGUS's output differs
between machines, so, as upstream, they check that every sequence is in the
alignment or tree. The guide-tree alignment uses a tree from
:func:`magus_guidetree`, as upstream's does.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.magus import magus_align, magus_guidetree
from tests.framework import env, fixture

SEQUENCES = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin-ref.fa"
SEATOXIN_IDS = ("1apf", "1ahl", "1atx", "1sh1", "1bds")


async def check_ids(file: File, name: str, marker: str = ">") -> None:
    assert file.name == name, f"output named {file.name}, expected {name}"
    async with file.open("rb") as fh:
        raw = bytes(await fh.read())
    text = (gzip.decompress(raw) if name.endswith(".gz") else raw).decode()
    missing = [i for i in SEATOXIN_IDS if f"{marker}{i}" not in text]
    assert not missing, f"{name} lacks {missing}"


@env.task
async def test_magus_align_seatoxin_uncompressed() -> None:
    out = await magus_align(await fixture(SEQUENCES), prefix="test")
    await check_ids(out, "test.aln")


@env.task
async def test_magus_align_seatoxin_compressed() -> None:
    out = await magus_align(await fixture(SEQUENCES), prefix="test", compress=True)
    await check_ids(out, "test.aln.gz")


@env.task
async def test_magus_align_seatoxin_guide_tree() -> None:
    sequences = await fixture(SEQUENCES)
    tree = await magus_guidetree(sequences, prefix="test")
    out = await magus_align(sequences, tree=tree, prefix="test", compress=True)
    await check_ids(out, "test.aln.gz")


@env.task
async def test_magus_guidetree_seatoxin() -> None:
    tree = await magus_guidetree(await fixture(SEQUENCES), prefix="test")
    await check_ids(tree, "test.tree", marker="")


tests = [
    test_magus_align_seatoxin_uncompressed,
    test_magus_align_seatoxin_compressed,
    test_magus_align_seatoxin_guide_tree,
    test_magus_guidetree_seatoxin,
]
