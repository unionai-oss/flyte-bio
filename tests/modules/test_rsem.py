"""Tests for flyte_bio.modules.rsem.

The case is upstream's rsem/preparereference test (modules commit 9a1cc9f,
RSEM 1.3.1), with its snapshot md5s for every file of the reference.
"""

from flyte_bio.modules.rsem import rsem_preparereference
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
TRANSCRIPTS_MD5 = "050c521a2719c2ae48267c1e65218f29"
INDEX_MD5S = {
    "genome.chrlist": "b190587cae0531f3cf25552d8aa674db",
    "genome.fasta": "f315020d899597c1b57e5fe9f60f4c3e",
    "genome.grp": "c2848a8b6d495956c11ec53efc1de67e",
    "genome.idx.fa": TRANSCRIPTS_MD5,
    "genome.n2g.idx.fa": TRANSCRIPTS_MD5,
    "genome.seq": "94da0c6b88c33e63c9a052a11f4f57c1",
    "genome.ti": "c9e4ae8d4d13a504eec2acf1b8589a66",
    "genome.transcripts.fa": TRANSCRIPTS_MD5,
}


@env.task
async def test_rsem_preparereference_homo_sapiens() -> None:
    # upstream case: rsem/preparereference "homo_sapiens"
    reference = await rsem_preparereference(
        await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.gtf")
    )
    files = {f.name: f for f in await reference.index.list_files()}
    assert sorted(files) == sorted(INDEX_MD5S), f"rsem reference holds {sorted(files)}"
    for name, md5 in INDEX_MD5S.items():
        await assert_md5(files[name], md5, label=name)
    assert reference.transcript_fasta.name == "genome.transcripts.fa", reference.transcript_fasta.name
    await assert_md5(reference.transcript_fasta, TRANSCRIPTS_MD5, label="transcript_fasta")


tests = [test_rsem_preparereference_homo_sapiens]
