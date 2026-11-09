"""Tests for flyte_bio.modules.ragtag.

The cases and expected md5s are upstream's ragtag/scaffold and ragtag/patch
tests (modules commit 9a1cc9f, RagTag 2.1.0), with the scaffold tests'
config args. nf-test computes the md5s of gzipped outputs decompressed. As
upstream, the alignment logs and the nucmer delta are checked by name.
"""

from flyte_bio.modules.ragtag import RagtagScaffoldResult, ragtag_patch, ragtag_scaffold
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta.gz"


async def scaffold(aligner: str) -> RagtagScaffoldResult:
    genome = await fixture(GENOME)
    return await ragtag_scaffold(genome, genome, prefix="test", args=f"-u --aligner {aligner}")


async def check_scaffolds(result: RagtagScaffoldResult) -> None:
    await assert_gunzipped_md5(result.corrected_assembly, "70b661fab5364a1c389972a771f97905", label="test.fasta.gz")
    await assert_md5(result.corrected_agp, "40fbf2d081c32880d8ce8187c529a80b", label="test.agp")
    await assert_md5(result.corrected_stats, "209e973e4bac1653b8d5fddb7fa13b63", label="test.stats")
    await assert_md5(result.confidence_txt, "4c7bc2c5128996cd72161c39fbfc3ce3", label="test.confidence.txt")


@env.task
async def test_ragtag_scaffold_minimap2() -> None:
    result = await scaffold("minimap2")
    await check_scaffolds(result)
    assert result.ragtag_paf is not None, "ragtag wrote no PAF"
    await assert_gunzipped_md5(result.ragtag_paf, "6bb30eea71e41b573e789e8d601dd778", label="test.asm.paf.gz")
    assert result.ragtag_paf_log is not None and result.ragtag_paf_log.name == "test.asm.paf.log"
    assert result.ragtag_delta is None and result.ragtag_delta_log is None


@env.task
async def test_ragtag_scaffold_nucmer() -> None:
    result = await scaffold("nucmer")
    await check_scaffolds(result)
    assert result.ragtag_paf is not None, "ragtag wrote no PAF"
    await assert_gunzipped_md5(result.ragtag_paf, "cb9db2cc6fef62f1f1aa3bc7700d94e6", label="test.asm.paf.gz")
    assert result.ragtag_delta is not None and result.ragtag_delta.name == "test.asm.delta.gz"
    assert result.ragtag_delta_log is not None and result.ragtag_delta_log.name == "test.asm.delta.log"
    assert result.ragtag_paf_log is None


@env.task
async def test_ragtag_patch() -> None:
    genome = await fixture(GENOME)
    result = await ragtag_patch(genome, genome, prefix="test")
    await assert_md5(result.patch_fasta, "5cf615df690061ab15e4fee62abf3ebc", label="test.patch.fasta")
    await assert_md5(result.patch_agp, "7878fd4e42ecb2bfccd7565d5ed6b625", label="test.patch.agp")
    await assert_md5(result.patch_components_fasta, "65bd2563dfc2564d5edf0e8d24257032", label="test.comps.fasta")
    await assert_md5(result.target_splits_agp, "ac3460a377daaf3e3ce37f499e561968", label="test.ctg.agp")
    await assert_md5(result.target_splits_fasta, "5cf615df690061ab15e4fee62abf3ebc", label="test.ctg.fasta")


tests = [test_ragtag_scaffold_minimap2, test_ragtag_scaffold_nucmer, test_ragtag_patch]
