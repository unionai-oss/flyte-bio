"""Tests for flyte_bio.modules.adapterremoval.

The cases and expected md5s are upstream's adapterremoval tests (modules
commit 9a1cc9f, AdapterRemoval 2.3.4): single-end, paired-end, paired-end
with ``--collapse``, and paired-end with upstream's adapter list. nf-test
hashes the gzipped reads decompressed. As upstream, the collapse case checks
that the (empty) quality-trimmed collapsed reads are written, by name.
"""

from flyte.io import File

from flyte_bio.modules.adapterremoval import AdapterRemovalResult, adapterremoval
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


async def adapterremoval_paired(args: str = "", adapter_list: File | None = None) -> AdapterRemovalResult:
    return await adapterremoval(
        await fixture(FASTQ + "test_1.fastq.gz"),
        reads_2=await fixture(FASTQ + "test_2.fastq.gz"),
        adapter_list=adapter_list,
        prefix="test",
        args=args,
    )


async def check_adapterremoval_pairs(result: AdapterRemovalResult, md5s: list[str]) -> None:
    names = [f.name for f in result.paired_truncated]
    assert names == ["test.pair1.truncated.fastq.gz", "test.pair2.truncated.fastq.gz"], names
    for file, md5 in zip(result.paired_truncated, md5s, strict=True):
        await assert_gunzipped_md5(file, md5, label=file.name)


@env.task
async def test_adapterremoval_single_end() -> None:
    # upstream case: adapterremoval "single-end - sarscov2 - [fastq]"
    result = await adapterremoval(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    assert result.singles_truncated is not None, "no test.truncated.fastq.gz"
    await assert_gunzipped_md5(result.singles_truncated, "119d1b1a0a71ca6e080ff7c53ee0b690", label="truncated")
    await assert_md5(result.settings, "9065cf46619a40116bdd55b697eaec33", label="test.settings")


@env.task
async def test_adapterremoval_paired_end() -> None:
    # upstream case: adapterremoval "paired-end - sarscov2 - [fastq]"
    result = await adapterremoval_paired()
    await check_adapterremoval_pairs(result, ["e3da014fbb9b428e952c62e8f0fb6402", "2ebae722295ea66d84075a3b042e2b42"])
    await assert_md5(result.settings, "9b7f6eeb6af14639588b76d434bd2d8f", label="test.settings")


@env.task
async def test_adapterremoval_paired_end_collapse() -> None:
    # upstream case: adapterremoval "paired-end collapse - sarscov2 - [fastq]"
    result = await adapterremoval_paired(args="--collapse")
    await check_adapterremoval_pairs(result, ["d6113ac35457dc942e4e47d6530e1d5e", "304c48e7ad50d46acf73ae6de4014f64"])
    assert result.collapsed is not None, "no test.collapsed.fastq.gz"
    await assert_gunzipped_md5(result.collapsed, "369452751050a7f1e31b839702d61417", label="collapsed")
    assert result.collapsed_truncated is not None, "no test.collapsed.truncated.fastq.gz"
    assert result.collapsed_truncated.name == "test.collapsed.truncated.fastq.gz", result.collapsed_truncated.name
    await assert_md5(result.settings, "6059588d1cac39c081d15c2024fcdf04", label="test.settings")


@env.task
async def test_adapterremoval_paired_end_adapterlist() -> None:
    # upstream case: adapterremoval "paired-end adapterlist - sarscov2 - [fastq]"
    adapters = await fixture("delete_me/adapterremoval/adapterremoval_adapterlist.txt")
    result = await adapterremoval_paired(adapter_list=adapters)
    await check_adapterremoval_pairs(result, ["e3da014fbb9b428e952c62e8f0fb6402", "2ebae722295ea66d84075a3b042e2b42"])
    await assert_md5(result.settings, "5eca9942f2e3a970737f81251c7f49e7", label="test.settings")


tests = [
    test_adapterremoval_single_end,
    test_adapterremoval_paired_end,
    test_adapterremoval_paired_end_collapse,
    test_adapterremoval_paired_end_adapterlist,
]
