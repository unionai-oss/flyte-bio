"""Tests for flyte_bio.modules.rmats.

The cases and expected md5s are upstream's rmats/prep tests (modules commit
9a1cc9f, rMATS 4.3.0), with their test config's ``ext.args``: every case
but "no clipping" (which loads no config) runs with
``--variable-read-length --allow-clipping``, plus ``-t single`` or the
strandedness's ``--libType``, and ``--novelSS --mil 48 --mel 488`` for the
novel-splice case. As upstream, the read-outcome counts are also checked:
reads are used, and none are clipped unless clipping is disallowed.
"""

import re

from flyte.io import File

from flyte_bio.modules.rmats import RmatsPrepResult, rmats_prep
from tests.framework import assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/test.rna.paired_end.bam"
GTF = "genomics/homo_sapiens/genome/genome.gtf"
CONFIG_ARGS = "--variable-read-length --allow-clipping"


async def prep(args: str) -> RmatsPrepResult:
    return await rmats_prep(await fixture(BAM), await fixture(GTF), read_length=150, prefix="test", args=args)


async def outcome(file: File, line: int, name: str) -> int:
    async with file.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    found = re.search(rf"{name}: (\d+)", lines[line])
    assert found, f"read outcomes line {line + 1} lacks {name}: {lines[line]!r}"
    return int(found.group(1))


async def check(result: RmatsPrepResult, rmats: str, read_outcomes: str) -> int:
    await assert_md5(result.rmats, rmats, label="rmats")
    await assert_md5(result.read_outcomes, read_outcomes, label="read outcomes")
    used = await outcome(result.read_outcomes, 1, "USED")
    assert used > 0, "no reads used"
    return used


@env.task
async def test_rmats_prep_paired_unstranded() -> None:
    # upstream case: rmats/prep "homo_sapiens - paired unstranded rmats prep"
    result = await prep(CONFIG_ARGS)
    await check(result, "3b847b81badea6cbf623d6072609a743", "18975bd7f2e4a3ef9186fc4f1ebca5f2")
    assert await outcome(result.read_outcomes, 9, "CLIPPED") == 0, "clipped reads despite --allow-clipping"


@env.task
async def test_rmats_prep_single_end_unstranded() -> None:
    # upstream case: rmats/prep "homo_sapiens - single-end unstranded rmats prep"
    result = await prep(f"{CONFIG_ARGS} -t single")
    await check(result, "3b847b81badea6cbf623d6072609a743", "18975bd7f2e4a3ef9186fc4f1ebca5f2")


@env.task
async def test_rmats_prep_paired_forward() -> None:
    # upstream case: rmats/prep "homo_sapiens - paired forward rmats prep"
    result = await prep(f"{CONFIG_ARGS} --libType fr-firststrand")
    await check(result, "042a54dec896ab6b00ddd13dd49218f5", "96bfd827a39f14566acde8adc1249145")


@env.task
async def test_rmats_prep_paired_reverse() -> None:
    # upstream case: rmats/prep "homo_sapiens - paired reverse rmats prep"
    result = await prep(f"{CONFIG_ARGS} --libType fr-secondstrand")
    await check(result, "a879b17464d0f1d6025c10db1148f881", "ebd1a7e14896503c2dce3d6ed6b866d3")


@env.task
async def test_rmats_prep_paired_novel_splice() -> None:
    # upstream case: rmats/prep "homo_sapiens - paired unstranded novel splice rmats prep"
    result = await prep(f"{CONFIG_ARGS} --novelSS --mil 48 --mel 488")
    await check(result, "3b847b81badea6cbf623d6072609a743", "18975bd7f2e4a3ef9186fc4f1ebca5f2")


@env.task
async def test_rmats_prep_paired_no_clipping() -> None:
    # upstream case: rmats/prep "homo_sapiens - paired unstranded no clipping rmats prep"
    result = await prep("")
    used = await check(result, "3800cad59b9d91d1da63d40dce08173f", "fb7d9dab371bb7097b9f2358e61170e0")
    clipped = await outcome(result.read_outcomes, 9, "CLIPPED")
    assert clipped > used, f"{clipped} reads clipped, expected more than the {used} used"


tests = [
    test_rmats_prep_paired_unstranded,
    test_rmats_prep_single_end_unstranded,
    test_rmats_prep_paired_forward,
    test_rmats_prep_paired_reverse,
    test_rmats_prep_paired_novel_splice,
    test_rmats_prep_paired_no_clipping,
]
