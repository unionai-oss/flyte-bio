"""Tests for flyte_bio.modules.pilon.

The cases and expected md5s are upstream's pilon tests (modules commit
9a1cc9f, Pilon 1.24): both modes write the same improved assembly and no
other output.
"""

from flyte_bio.modules.pilon import PilonResult, pilon
from tests.framework import assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"


async def run_pilon(mode: str) -> PilonResult:
    fasta = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    return await pilon(fasta, await fixture(BAM), await fixture(BAM + ".bai"), pilon_mode=mode, prefix="test")


async def check(result: PilonResult) -> None:
    assert result.improved_assembly.name == "test.fasta", f"Pilon wrote {result.improved_assembly.name}"
    await assert_md5(result.improved_assembly, "2e881994820a5a641da9ea594ab4958f", label="test.fasta")
    assert result.vcf is None and result.change_record is None and result.tracks == [], "unexpected Pilon outputs"


@env.task
async def test_pilon_bam() -> None:
    await check(await run_pilon("bam"))


@env.task
async def test_pilon_frags() -> None:
    await check(await run_pilon("frags"))


tests = [test_pilon_bam, test_pilon_frags]
