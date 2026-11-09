"""Tests for flyte_bio.modules.cooler.

The cases are upstream's cooler tests (modules commit 9a1cc9f, cooler 0.10.4;
dump 0.10.3). The digest and makebins BEDs and the dump of a toy cooler are
pure data, so they're checked against upstream's md5s. ``.cool`` and ``.mcool`` files are
HDF5 and record their creation date, so upstream snapshots only their names;
these tests check the names too, then that the file is HDF5 and, through
:func:`cooler_dump`, that it holds what the step should have written.
"""

from flyte.io import File

from flyte_bio.modules.cooler import (
    CoolerResult,
    cooler_balance,
    cooler_cload,
    cooler_digest,
    cooler_dump,
    cooler_makebins,
    cooler_zoomify,
    file_name,
)
from tests.framework import assert_md5, env, fixture

CLOAD = "genomics/homo_sapiens/cooler/cload/hg19/"
TOY = "genomics/homo_sapiens/cooler/merge/toy/toy.symm.upper.2.cool"
DUMP_TOY = "https://raw.githubusercontent.com/open2c/cooler/master/tests/data/toy.asymm.16.cool"
HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check_cooler(result: CoolerResult, name: str) -> None:
    assert file_name(result.output) == name, f"cooler wrote {file_name(result.output)}, expected {name}"
    async with result.output.open("rb") as fh:
        head = bytes(await fh.read(len(HDF5_SIGNATURE)))
    assert head == HDF5_SIGNATURE, f"{name} is not HDF5"


async def check_pixels(cool: File, label: str, resolution: int | None = None) -> None:
    dumped = await cooler_dump(cool, resolution=resolution, prefix="pixels")
    lines = (await text(dumped.output)).splitlines()
    assert lines, f"{label}: no pixels"
    assert all(len(line.split("\t")) == 3 for line in lines), f"{label}: pixels aren't bin1, bin2, count"


@env.task
async def test_cooler_cload_pairix() -> None:
    result = await cooler_cload(
        await fixture(CLOAD + "hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz"),
        await fixture(CLOAD + "hg19.chrom.sizes"),
        mode="pairix",
        bin_size=2000000,
        index=await fixture(CLOAD + "hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz.px2"),
        prefix="test_pairix",
    )
    await check_cooler(result, "test_pairix.cool")
    await check_pixels(result.output, "pairix")


@env.task
async def test_cooler_cload_pairs() -> None:
    result = await cooler_cload(
        await fixture(CLOAD + "hg19.sample1.pairs"),
        await fixture(CLOAD + "hg19.chrom.sizes"),
        mode="pairs",
        bin_size=2000000,
        prefix="test_pairs",
        args="--chrom1 1 --pos1 2 --chrom2 4 --pos2 5 -N",
    )
    await check_cooler(result, "test_pairs.cool")
    await check_pixels(result.output, "pairs")


@env.task
async def test_cooler_cload_tabix() -> None:
    result = await cooler_cload(
        await fixture(CLOAD + "hg19.GM12878-MboI.pairs.subsample.sorted.possrt.txt.gz"),
        await fixture(CLOAD + "hg19.chrom.sizes"),
        mode="tabix",
        bin_size=2000000,
        index=await fixture(CLOAD + "hg19.GM12878-MboI.pairs.subsample.sorted.possrt.txt.gz.tbi"),
        prefix="test_tabix",
    )
    await check_cooler(result, "test_tabix.cool")
    await check_pixels(result.output, "tabix")


@env.task
async def test_cooler_digest() -> None:
    result = await cooler_digest(
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        await fixture("genomics/sarscov2/genome/genome.sizes"),
        "CviQI",
    )
    assert file_name(result.output) == "genome_CviQI.bed", file_name(result.output)
    await assert_md5(result.output, "f02d3cfaf642b27d4cc197209845e096", label="digest")


@env.task
async def test_cooler_makebins() -> None:
    result = await cooler_makebins(await fixture(CLOAD + "hg19.chrom.sizes"), 1000000, prefix="test")
    assert file_name(result.output) == "test.bed"
    await assert_md5(result.output, "99e887ba1237ebcdcf31954476ab0adf", label="makebins")


@env.task
async def test_cooler_balance() -> None:
    result = await cooler_balance(await fixture(TOY), prefix="test")
    await check_cooler(result, "test.cool")
    bins = await cooler_dump(result.output, prefix="bins", args="-t bins --header")
    header = (await text(bins.output)).splitlines()[0].split("\t")
    assert "weight" in header, f"balanced bins lack a weight column: {header}"


@env.task
async def test_cooler_dump() -> None:
    result = await cooler_dump(await fixture(DUMP_TOY), prefix="test")
    assert file_name(result.output) == "test.bedpe"
    await assert_md5(result.output, "38e9b0b8cc74f55a15e8ab01023048bd", label="dump")


@env.task
async def test_cooler_zoomify() -> None:
    result = await cooler_zoomify(await fixture(TOY), prefix="test", args="-r 2,4,8")
    await check_cooler(result, "test.mcool")
    for resolution in (2, 4, 8):
        await check_pixels(result.output, f"zoomify {resolution}", resolution=resolution)


tests = [
    test_cooler_cload_pairix,
    test_cooler_cload_pairs,
    test_cooler_cload_tabix,
    test_cooler_digest,
    test_cooler_makebins,
    test_cooler_balance,
    test_cooler_dump,
    test_cooler_zoomify,
]
