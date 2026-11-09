"""Tests for flyte_bio.modules.liftoff.

The case and expected values are upstream's liftoff test (modules commit
9a1cc9f, Liftoff 1.6.3) with its test config's ``-polish``: chr1's GTF lifted
from chr1 onto chr21. As upstream, the unmapped features are compared by md5
and both annotations must contain a known lifted exon.
"""

from flyte.io import File

from flyte_bio.modules.liftoff import liftoff
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
EXON = "chr21\tLiftoff\texon\t34608061\t34608118\t.\t+\t."


async def liftoff_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_liftoff_chr1_to_chr21() -> None:
    # upstream case: "homo_sapiens-genome_21_fasta-genome_1_fasta-genome_1_gtf"
    result = await liftoff(
        await fixture(GENOME + "chr21/sequence/genome.fasta"),
        await fixture(GENOME + "chr1/genome.fasta.gz"),
        ref_annotation=await fixture(GENOME + "chr1/genome.gtf"),
        prefix="test",
        args="-polish",
    )
    await assert_md5(result.unmapped_txt, "7391d10df6e15db356b084c9af5259e4", label="unmapped.txt")
    assert EXON in await liftoff_text(result.gff3), "gff3 lacks the lifted exon"
    assert result.polished_gff3 is not None, "liftoff wrote no polished gff3"
    assert EXON in await liftoff_text(result.polished_gff3), "polished gff3 lacks the lifted exon"


tests = [test_liftoff_chr1_to_chr21]
