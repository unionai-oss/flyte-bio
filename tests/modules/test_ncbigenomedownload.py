"""Tests for flyte_bio.modules.ncbigenomedownload.

The case and expected values are upstream's ncbigenomedownload test (modules
commit 9a1cc9f, ncbi-genome-download 0.3.3): GenBank, FASTA and assembly
stats of one Staphylococcus aureus assembly, downloaded from NCBI at run
time. As upstream's snapshot: lines 4-8 of the GenBank and FASTA files.

Upstream also snapshots the assembly stats' md5 (23b5b3a7…, March 2026), but
NCBI has since rewritten that file (Last-Modified 22 May 2026): what the tool
downloads now is byte-identical to NCBI's current file (f9f9aec3…), not to
the snapshot. So the stats are checked by content: the assembly's name and
accession, and its total length.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.ncbigenomedownload import ncbigenomedownload
from tests.framework import env

GBK_LINES = [
    "ACCESSION   NC_007795",
    "VERSION     NC_007795.1",
    "DBLINK      BioProject: PRJNA57795",
    "            Assembly: GCF_000013425.1",
    "KEYWORDS    RefSeq.",
]
FNA_LINES = [
    "GATACTCTATAGAAGGAAAAGTTATCCACTTATGCACATTTATAGTTTTCAGAATTGTGGATAATTAGAAATTACACACA",
    "AAGTTATACTATTTTTAGCAACATATTCACAGGTATTTGACATATAGAGAACTGAAAAAGTATAATTGTGTGGATAAGTC",
    "GTCCAACTCATGATTTTATAAGGATTTATTTATTGATTTTTACATAAAAATACTGTGCATAACTAATAAGCAAGATAAAG",
    "TTATCCACCGATTGTTATTAACTTGTGGATAATTATTAACATGGTGTGTTTAGAAGTTATCCACGGCTGTTATTTTTGTG",
    "TATAACTTAAAAATTTAAGAAAGATGGAGTAAATTTATGTCGGAAAAAGAAATTTGGGAAAAAGTGCTTGAAATTGCTCA",
]

STATS_LINES = [
    "# Assembly name:  ASM1342v1",
    "# RefSeq assembly accession: GCF_000013425.1",
    "all\tall\tall\tall\ttotal-length\t2821361",
]


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def gunzipped_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode().splitlines()


@env.task
async def test_ncbigenomedownload_bacteria() -> None:
    result = await ncbigenomedownload("bacteria", args="-A GCF_000013425.1 --formats genbank,fasta,assembly-stats -N")
    assert len(result.gbk) == 1 and len(result.fna) == 1, (result.gbk, result.fna)
    assert (await gunzipped_lines(result.gbk[0]))[3:8] == GBK_LINES, "GenBank lines 4-8 differ"
    assert (await gunzipped_lines(result.fna[0]))[3:8] == FNA_LINES, "FASTA lines 4-8 differ"
    assert [f.name for f in result.stats] == ["GCF_000013425.1_ASM1342v1_assembly_stats.txt"], result.stats
    stats = await text(result.stats[0])
    for line in STATS_LINES:
        assert line in stats.splitlines(), f"assembly stats lack {line!r}"


tests = [test_ncbigenomedownload_bacteria]
