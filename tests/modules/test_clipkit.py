"""Tests for flyte_bio.modules.clipkit.

The cases and expected md5s are upstream's clipkit tests (modules commit
9a1cc9f, ClipKIT 2.11.4), on the PF14720 seed alignment decompressed as
upstream's gunzip setup does. Upstream snapshots the alignment and log by
md5, and the console output and JSON report by line count.
"""

import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.clipkit import ClipkitResult, clipkit
from tests.framework import assert_md5, env, fixture

SEED = "delete_me/hmmer/PF14720_seed.alnfaa.gz"


async def file_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


def local_file(name: str, data: bytes) -> str:
    path = Path(tempfile.mkdtemp()) / name
    path.write_bytes(data)
    return str(path)


async def seed_alignment() -> File:
    async with (await fixture(SEED)).open("rb") as fh:
        data = gzip.decompress(bytes(await fh.read()))
    return await File.from_local(local_file("PF14720_seed.alnfaa", data))


async def check(result: ClipkitResult, name: str, alignment: str, log: str) -> None:
    assert result.alignment.name == name, result.alignment.name
    await assert_md5(result.alignment, alignment, label=name)
    await assert_md5(result.log, log, label=f"{name}.log")
    lines = len((await file_text(result.metadata)).splitlines())
    assert lines == 34, f"{name}.txt has {lines} lines, expected 34"


@env.task
async def test_clipkit_aln() -> None:
    result = await clipkit(await seed_alignment(), prefix="PF14720")
    await check(result, "PF14720.clipkit", "97e4328b5990f5e5a6241c225ed6d922", "53bb80178adeb0fa6b4df882a6963c28")
    assert result.complementary is None and result.json is None


@env.task
async def test_clipkit_gap_threshold_ends_only() -> None:
    result = await clipkit(
        await seed_alignment(), out_format="alnfaa", prefix="PF14720", args="--mode gappy --gaps 0.1 --ends_only"
    )
    await check(result, "PF14720.alnfaa", "c1ab1d4f55e928c24d03dd1fba1d8f27", "ead8293153a057d4e564e45838fb39e0")


@env.task
async def test_clipkit_complementary_json() -> None:
    result = await clipkit(await seed_alignment(), prefix="PF14720", args="--complementary --report_json")
    await check(result, "PF14720.clipkit", "97e4328b5990f5e5a6241c225ed6d922", "53bb80178adeb0fa6b4df882a6963c28")
    assert result.complementary is not None and result.json is not None
    await assert_md5(result.complementary, "3fda46ab66fb05a3ec36b7eda4c285c0", label="PF14720.clipkit.complement")
    lines = len((await file_text(result.json)).splitlines())
    assert lines == 31, f"PF14720.clipkit.report.json has {lines} lines, expected 31"


@env.task
async def test_clipkit_auxiliary_file_cst() -> None:
    sites = "".join(f"{site}\n" for site in ["1\tkeep", "2\ttrim", "3\tkeep", "4\ttrim", "5\tkeep"])
    aux = await File.from_local(local_file("cst_auxiliary.tsv", sites.encode()))
    result = await clipkit(await seed_alignment(), auxiliary_file=aux, prefix="PF14720", args="--mode cst")
    await check(result, "PF14720.clipkit", "cec05fd94f4063cd5d4dd3630d78b0c9", "67ed60637c2adfa059a34f413418d6b3")


tests = [
    test_clipkit_aln,
    test_clipkit_gap_threshold_ends_only,
    test_clipkit_complementary_json,
    test_clipkit_auxiliary_file_cst,
]
