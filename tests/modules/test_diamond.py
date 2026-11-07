"""Tests for flyte_bio.modules.diamond.

The cases are upstream's diamond/cluster tests (modules commit 9a1cc9f,
DIAMOND 2.2.2): the human reference proteome, as a FASTA and as a ``.dmnd``
database, with no ``args``. As upstream, they check the cluster table's line
count (one line per protein), not its md5.
"""

from flyte.io import File

from flyte_bio.modules.diamond import diamond_cluster
from tests.framework import env, fixture

PROTEOME = "proteomics/database/UP000005640_9606"


async def line_count(file: File) -> int:
    async with file.open("rb") as fh:
        return len(bytes(await fh.read()).decode().splitlines())


async def check(db: File) -> None:
    table = await diamond_cluster(db, prefix="human")
    assert table.name == "human.tsv", f"cluster table is {table.name}"
    lines = await line_count(table)
    assert lines == 20610, f"human.tsv has {lines} lines, expected 20610"


@env.task
async def test_diamond_cluster_human_fasta() -> None:
    await check(await fixture(f"{PROTEOME}.fasta"))


@env.task
async def test_diamond_cluster_human_dmnd() -> None:
    await check(await fixture(f"{PROTEOME}.dmnd"))


tests = [
    test_diamond_cluster_human_fasta,
    test_diamond_cluster_human_dmnd,
]
