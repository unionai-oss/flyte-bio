"""Tests for flyte_bio.modules.porechop.

The cases and expected md5s are upstream's porechop/porechop (Porechop 0.2.4)
and porechop/abi (Porechop_ABI 0.5.0) tests, modules commit 9a1cc9f. As
upstream's porechop test config, Porechop's prefix is ``test_porechop``;
Porechop_ABI keeps its default ``test.porechop_abi``. nf-test computes the
md5s over the decompressed reads; the logs are unstable, so only their
presence is checked.
"""

import os

from flyte.io import File

from flyte_bio.modules.porechop import PorechopResult, porechop, porechop_abi
from tests.framework import assert_gunzipped_md5, env, fixture

NANOPORE = "genomics/sarscov2/nanopore/fastq/test.fastq.gz"
TRIMMED_MD5 = "886fdb859fb50e0dddd35007bcff043e"


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def check(result: PorechopResult, prefix: str) -> None:
    assert file_name(result.reads) == f"{prefix}.fastq.gz", f"reads named {file_name(result.reads)}"
    assert file_name(result.log) == f"{prefix}.log", f"log named {file_name(result.log)}"
    await assert_gunzipped_md5(result.reads, TRIMMED_MD5, label="reads")


@env.task
async def test_porechop_porechop_nanopore() -> None:
    await check(await porechop(await fixture(NANOPORE), prefix="test_porechop"), "test_porechop")


@env.task
async def test_porechop_abi_nanopore() -> None:
    await check(await porechop_abi(await fixture(NANOPORE), prefix="test.porechop_abi"), "test.porechop_abi")


tests = [test_porechop_porechop_nanopore, test_porechop_abi_nanopore]
