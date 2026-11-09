"""Tests for flyte_bio.modules.links.

The cases and expected md5s are upstream's links tests (modules commit
9a1cc9f, LINKS 2.0.1), with their test config's distances. As upstream, the
log and scaffold graph are checked by name, the pairing issues to be empty
and, scaffolding scaffolds, the tig-pair checkpoint to be empty.
"""

from flyte.io import File

from flyte_bio.modules.links import LinksResult, links
from tests.framework import assert_md5, env, fixture

ARGS = "-d 1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000"
FASTA = "genomics/sarscov2/illumina/fasta/"
BLOOM = "23737e49d9a2f070b312da844201b494"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: LinksResult, scaffolds: str, scaffolds_fa: str, correspondence: str) -> None:
    assert result.log.name == "test.log" and result.scaffolds_graph.name == "test.gv"
    assert result.tigpair_checkpoint.name == "test.tigpair_checkpoint.tsv"
    assert await text(result.pairing_issues) == "", "LINKS reported pairing issues"
    await assert_md5(result.scaffolds_csv, scaffolds, label="test.scaffolds")
    await assert_md5(result.scaffolds_fasta, scaffolds_fa, label="test.scaffolds.fa")
    await assert_md5(result.bloom, BLOOM, label="test.bloom")
    await assert_md5(result.assembly_correspondence, correspondence, label="test.assembly_correspondence.tsv")


@env.task
async def test_links_contigs() -> None:
    contigs = await fixture(FASTA + "contigs.fasta")
    result = await links(contigs, [await fixture("genomics/sarscov2/genome/genome.fasta")], prefix="test", args=ARGS)
    await check(
        result,
        scaffolds="41c129edd1e66140fcfb7efce81197ad",
        scaffolds_fa="8abc4f609d0ad415f900b0046b38a72b",
        correspondence="b36e951b0a1bb4b1c1ccd50925392e3d",
    )


@env.task
async def test_links_scaffolds() -> None:
    scaffolds = await fixture(FASTA + "scaffolds.fasta")
    result = await links(scaffolds, [await fixture(FASTA + "contigs.fasta")], prefix="test", args=ARGS)
    await check(
        result,
        scaffolds="095cc323b3af3a7873c8b80cf3736a1f",
        scaffolds_fa="b8c7938abbc3d2f9b5c3d709d43b4a60",
        correspondence="a65d30663dce705d382df52ab87ca8a4",
    )
    assert await text(result.tigpair_checkpoint) == "", "LINKS wrote tig pairs"


tests = [test_links_contigs, test_links_scaffolds]
