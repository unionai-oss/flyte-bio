"""Tests for flyte_bio.modules.yahs.

The cases and expected md5s are upstream's yahs tests (modules commit
9a1cc9f, YaHS 1.2.2). The second rescaffolds from the first's final AGP with
``--no-contig-ec``.
"""

from flyte.io import File

from flyte_bio.modules.yahs import YahsResult, yahs
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/genome.fasta"
BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"


async def inputs() -> tuple[File, File, File]:
    return await fixture(GENOME), await fixture(GENOME + ".fai"), await fixture(BAM)


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: YahsResult, initial_break: bool) -> None:
    assert result.scaffolds_fasta is not None and result.scaffolds_agp is not None, "YaHS wrote no final scaffolds"
    await assert_md5(result.scaffolds_fasta, "a767604036a4cd6980ebd24f11e4dd95", label="test_scaffolds_final.fa")
    await assert_md5(result.scaffolds_agp, "374235b079d9e1c738c6f12697029b78", label="test_scaffolds_final.agp")
    names = [f.name for f in result.initial_break_agp]
    if initial_break:
        assert names == ["test_inital_break_01.agp"], f"initial-break AGPs {names}"
        await assert_md5(result.initial_break_agp[0], "374235b079d9e1c738c6f12697029b78", label=names[0])
    else:
        assert names == [], f"unexpected initial-break AGPs {names}"
    assert result.round_agp == [], f"unexpected round AGPs {[f.name for f in result.round_agp]}"
    await assert_md5(result.binary, "da45a2ec8a97fc24783e9a63373db379", label="test.bin")
    assert "Real time" in (await text(result.log)).splitlines()[-1], "YaHS log doesn't end with its run time"


@env.task
async def test_yahs_bam() -> None:
    await check(await yahs(*await inputs(), prefix="test"), initial_break=True)


@env.task
async def test_yahs_bam_agp() -> None:
    initial = await yahs(*await inputs(), prefix="test")
    assert initial.scaffolds_agp is not None, "YaHS wrote no final AGP"
    result = await yahs(*await inputs(), agp=initial.scaffolds_agp, prefix="test", args="--no-contig-ec")
    await check(result, initial_break=False)


tests = [test_yahs_bam, test_yahs_bam_agp]
