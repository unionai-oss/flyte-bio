"""Tests for flyte_bio.modules.pairtools.

The cases are upstream's pairtools parse, sort, dedup, flip, restrict and
select tests (modules commit 9a1cc9f, pairtools 1.1.3), with their test
configs' prefixes and select condition. Each upstream case starts from its own
fixture rather than chaining the steps, and these do the same. pairtools
records its command line in the output headers, so upstream snapshots lines
3 to 7 of each decompressed output (header lines before the ``@PG`` records)
and the md5 of the statistics files; these tests check the same.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.pairtools import (
    pairtools_dedup,
    pairtools_flip,
    pairtools_parse,
    pairtools_restrict,
    pairtools_select,
    pairtools_sort,
)
from tests.framework import assert_md5, env, fixture

PAIRTOOLS = "genomics/homo_sapiens/pairtools/"
SELECT_CONDITION = "(pair_type == 'RU') or (pair_type == 'UR') or (pair_type == 'UU')"
FLIP_LINES = [
    "#chromosomes: chr1 chr2",
    "#chromsize: chr1 10000",
    "#chromsize: chr2 10000",
    "#samheader: @SQ\tSN:chr1\tLN:10000",
    "#samheader: @SQ\tSN:chr2\tLN:10000",
]
SELECT_LINES = [
    "#chromosomes: chr2 chr3 chr1",
    "#chromsize: chr2 100",
    "#chromsize: chr3 100",
    "#chromsize: chr1 100",
    "#samheader: @SQ\tSN:chr1\tLN:100",
]


async def check_lines(file: File, name: str, expected: list[str]) -> None:
    """``file`` is named ``name`` and its decompressed lines 3 to 7 are ``expected``."""
    assert file.name == name, f"pairtools wrote {file.name}, expected {name}"
    async with file.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    assert lines[3:8] == expected, f"{name}: lines 3-7 are {lines[3:8]}, expected {expected}"


@env.task
async def test_pairtools_parse() -> None:
    result = await pairtools_parse(
        await fixture(PAIRTOOLS + "mock.sam"), await fixture(PAIRTOOLS + "mock.chrom.sizes"), prefix="test.raw"
    )
    await check_lines(
        result.pairs,
        "test.raw.pairsam.gz",
        [
            "#chromsize: chr1 1000",
            "#chromsize: chr2 1000",
            "#samheader: @SQ\tSN:chr1\tLN:1000",
            "#samheader: @SQ\tSN:chr2\tLN:1000",
            "#samheader: @PG\tID:mock\tPN:mock\tVN:0.0.0\tCL:mock",
        ],
    )
    assert result.stat.name == "test.raw.pairsam.stat", result.stat.name
    await assert_md5(result.stat, "bdca4672c0f0325d6882b8417651334b", label="parse stat")


@env.task
async def test_pairtools_sort() -> None:
    sorted_pairs = await pairtools_sort(await fixture(PAIRTOOLS + "mock.pairsam"), prefix="test.sorted")
    await check_lines(
        sorted_pairs,
        "test.sorted.pairs.gz",
        [
            "#genome_assembly: unknown",
            "#chromosomes: : chr1 chr2 chr3",
            "#chromsize: chr2 100",
            "#chromsize: chr3 100",
            "#chromsize: chr1 100",
        ],
    )


@env.task
async def test_pairtools_dedup() -> None:
    result = await pairtools_dedup(await fixture(PAIRTOOLS + "mock.4dedup.pairsam"), prefix="test.dedup")
    await check_lines(
        result.pairs,
        "test.dedup.pairs.gz",
        [
            "#genome_assembly: unknown",
            "#chromosomes: chr1 chr2",
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type sam1 sam2",
            "readid3\tchr1\t1\tchr1\t20\t+\t+\tUU\t.\t.",
            "readid5\tchr1\t1\tchr1\t25\t+\t+\tUU\t.\t.",
        ],
    )
    assert result.stat.name == "test.dedup.pairs.stat", result.stat.name
    await assert_md5(result.stat, "d31a54a0a28e2687c6ece6ddf970dd5a", label="dedup stat")


@env.task
async def test_pairtools_flip() -> None:
    flipped = await pairtools_flip(
        await fixture(PAIRTOOLS + "mock.4flip.pairs"), await fixture(PAIRTOOLS + "mock.chrom.sizes"), prefix="test"
    )
    await check_lines(flipped, "test.flip.gz", FLIP_LINES)


@env.task
async def test_pairtools_restrict() -> None:
    restricted = await pairtools_restrict(
        await fixture(PAIRTOOLS + "mock.4flip.pairs"), await fixture(PAIRTOOLS + "frag.bed"), prefix="test.restrict"
    )
    await check_lines(restricted, "test.restrict.pairs.gz", FLIP_LINES)


@env.task
async def test_pairtools_select() -> None:
    result = await pairtools_select(await fixture(PAIRTOOLS + "mock.pairsam"), SELECT_CONDITION, prefix="test")
    await check_lines(result.selected, "test.selected.pairs.gz", SELECT_LINES)
    await check_lines(result.unselected, "test.unselected.pairs.gz", SELECT_LINES)


tests = [
    test_pairtools_parse,
    test_pairtools_sort,
    test_pairtools_dedup,
    test_pairtools_flip,
    test_pairtools_restrict,
    test_pairtools_select,
]
