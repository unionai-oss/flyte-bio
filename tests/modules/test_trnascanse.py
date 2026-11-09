"""Tests for flyte_bio.modules.trnascanse.

The case and expected md5s are upstream's tRNAscan-SE test (modules commit
9a1cc9f, tRNAscan-SE 2.0.13): every output written for the Bacteroides
fragilis genome. As upstream, the log and stats (timestamped) are checked by
content; the module's 4 CPUs match upstream's test CPUs.

Upstream's ``fasta`` output is always empty (it gzips the FASTA, then looks
for ``*.fasta``), so its snapshot has no FASTA md5; here the gzipped FASTA is
returned and checked to be nonempty.
"""

from flyte.io import File

from flyte_bio.modules.trnascanse import trnascanse
from tests.framework import assert_md5, assert_nonempty, env, fixture


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_trnascanse_write_all() -> None:
    result = await trnascanse(
        await fixture("genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz"),
        write_stats=True,
        write_fasta=True,
        write_gff=True,
        write_bed=True,
        prefix="test",
    )
    await assert_md5(result.tsv, "23c19ff33ce36f041803b5b3e7b7681b", label="tsv")
    assert result.gff is not None and result.bed is not None and result.stats is not None
    await assert_md5(result.gff, "7e91ca1cbf554e63fa2d31e78e4fe7dd", label="gff")
    await assert_md5(result.bed, "4ab4ec79d400145086474745569dc8ac", label="bed")
    first = (await text(result.log)).splitlines()[0]
    assert "Application: tRNAscan-SE" in first, f"log starts {first!r}"
    assert "scan results" in await text(result.stats), "stats lack 'scan results'"
    assert result.fasta is not None, "tRNAscan-SE wrote no FASTA"
    await assert_nonempty(result.fasta, label="fasta")


tests = [test_trnascanse_write_all]
