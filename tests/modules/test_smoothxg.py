"""Tests for flyte_bio.modules.smoothxg.

The cases are upstream's smoothxg tests (modules commit 9a1cc9f, smoothxg
0.8.0) with their test config's ``ext.args``,
``--write-msa-in-maf-format test.maf --n-haps 10``. Upstream snapshots only
the tool version and asserts a smoothed graph was written; these tests also
check the graph is a GFA with paths and that the MAF was written.
"""

from flyte_bio.modules.smoothxg import SmoothxgResult, smoothxg
from tests.framework import env, fixture

ARGS = "--write-msa-in-maf-format test.maf --n-haps 10"


async def check_smoothed(result: SmoothxgResult) -> None:
    assert result.gfa.path.endswith("/test.smoothxg.gfa"), result.gfa.path
    async with result.gfa.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert lines and lines[0].startswith("H\t"), lines[:1]
    assert any(line.startswith("S\t") for line in lines), "smoothed graph has no segments"
    assert any(line.startswith("P\t") for line in lines), "smoothed graph has no paths"
    assert result.maf is not None and result.maf.path.endswith("/test.maf"), result.maf


@env.task
async def test_smoothxg_assembly() -> None:
    # upstream case: smoothxg "sarscov2 - illumina - assembly_gfa"
    gfa = await fixture("genomics/sarscov2/illumina/gfa/assembly.gfa")
    await check_smoothed(await smoothxg(gfa, prefix="test", args=ARGS))


@env.task
async def test_smoothxg_pangenome() -> None:
    # upstream case: smoothxg "homo_sapiens - pangenome - pangenome_seqwish_gfa"
    gfa = await fixture("pangenomics/homo_sapiens/pangenome.seqwish.gfa")
    await check_smoothed(await smoothxg(gfa, prefix="test", args=ARGS))


tests = [test_smoothxg_assembly, test_smoothxg_pangenome]
