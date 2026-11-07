"""Tests for flyte_bio.modules.amulety.

The cases and expected md5s are upstream's amulety/translate and amulety/embed
tests (modules commit 9a1cc9f, amulety 2.1.3): translate a 10x PBMC BCR AIRR
rearrangement table with the test IgBLAST reference, then embed its heavy
chains with AntiBERTy. Upstream extracts the reference with untar first; here
the archive is passed to the task as it is.

The embedding's metadata table matches upstream's md5, but the embeddings
themselves don't: their float32 values differ in the last digits from
upstream's (the same here on 1 CPU, as upstream's test config sets, and on 4,
so it's the hardware's arithmetic, not threading). So the embeddings table is
checked by content: its header, one 512-dimension row per heavy chain of the
metadata, in the same order, and finite values.
"""

import math

from flyte.io import File

from flyte_bio.modules.amulety import amulety_embed, amulety_translate
from tests.framework import assert_md5, env, fixture

AIRR = "genomics/homo_sapiens/10xgenomics/cellranger/sc5p_v2_hs_PBMC_1k_bcr/sc5p_v2_hs_PBMC_1k_b_airr_rearrangement.tsv"
IGBLAST = "genomics/homo_sapiens/genome/index/igblast/igblast_base.tar"
ANTIBERTY_DIMENSIONS = 512
HEAVY_CHAINS = 82


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def translated() -> File:
    return await amulety_translate(await fixture(AIRR), await fixture(IGBLAST), prefix="test")


@env.task
async def test_amulety_translate() -> None:
    await assert_md5(await translated(), "70447131da302612108bdf10fc5cecef", label="test_translated.tsv")


@env.task
async def test_amulety_embed() -> None:
    result = await amulety_embed(await translated(), chain="H", model="antiberty", prefix="test")
    await assert_md5(result.metadata, "233805ee6b8c5ba151b1d864b2d02f4f", label="test_metadata.tsv")

    metadata = [line.split("\t") for line in (await text(result.metadata)).splitlines()]
    sequence_ids = [row[metadata[0].index("sequence_id")] for row in metadata[1:]]
    rows = [line.split("\t") for line in (await text(result.embedding)).splitlines()]
    header = ["sequence_id"] + [f"dim_{i}" for i in range(1, ANTIBERTY_DIMENSIONS + 1)]
    assert rows[0] == header, f"test.tsv header: {rows[0][:3]}... ({len(rows[0])} columns)"
    assert len(rows) - 1 == HEAVY_CHAINS, f"test.tsv has {len(rows) - 1} rows, expected {HEAVY_CHAINS}"
    assert [row[0] for row in rows[1:]] == sequence_ids, "test.tsv's sequences differ from its metadata's"
    for row in rows[1:]:
        assert len(row) == len(header), f"{row[0]}: {len(row)} columns"
        assert all(math.isfinite(float(v)) for v in row[1:]), f"{row[0]}: non-finite embedding"


tests = [test_amulety_translate, test_amulety_embed]
