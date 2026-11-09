"""Tests for flyte_bio.modules.wfmash.

The case is upstream's wfmash test (modules commit 9a1cc9f, wfmash 0.13.0)
with its test config's ``ext.args``, ``-s 500 -X``: an all-to-all mapping of
the bgzipped pangenome FASTA. Upstream marks the PAF unstable and checks only
its name; this test also checks it holds PAF records between the input's
sequences.
"""

from flyte.io import File

from flyte_bio.modules.wfmash import wfmash
from tests.framework import env, fixture

PANGENOME = "pangenomics/homo_sapiens/"


async def wfmash_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_wfmash_pangenome() -> None:
    # upstream case: wfmash "homo_sapiens - pangenome - pangenome_fa_bgzip"
    fasta = await fixture(PANGENOME + "pangenome.fa.gz")
    gzi = await fixture(PANGENOME + "pangenome.fa.gz.gzi")
    fai = await fixture(PANGENOME + "pangenome.fa.gz.fai")
    result = await wfmash(fasta, gzi, fai, query_self=True, prefix="test", args="-s 500 -X")
    assert result.paf.path.endswith("/test.paf"), result.paf.path
    names = {line.split("\t")[0] for line in (await wfmash_text(fai)).splitlines()}
    records = [line.split("\t") for line in (await wfmash_text(result.paf)).splitlines()]
    assert records, "wfmash wrote no mappings"
    for fields in records:
        assert len(fields) >= 12 and fields[0] in names and fields[5] in names, fields[:6]


tests = [test_wfmash_pangenome]
