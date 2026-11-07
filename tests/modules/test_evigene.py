"""Tests for flyte_bio.modules.evigene.

The case and expected md5s are upstream's evigene/tr2aacds test (modules
commit 9a1cc9f, EvidentialGene 23.7.15). As upstream's snapshot: the dropset
files' md5s, and the okayset's file names. The module's 4 CPUs match
upstream's test CPUs.
"""

from flyte_bio.modules.evigene import evigene_tr2aacds
from tests.framework import assert_md5, env, fixture

DROPSET = {
    "test.drop.aa": "b43c660c34a419ec35090c163a4e4a54",
    "test.drop.cds": "5adb810e379fde7d2057a34ead7eaa58",
    "test.drop.tr": "5f38aa1a7e7adb09ce34c0bd74727d67",
}
OKAYSET = [
    "test.ann.txt",
    "test.cull.aa",
    "test.cull.cds",
    "test.cull.mrna",
    "test.genesum.txt",
    "test.mainalt.tab",
    "test.okay.aa",
    "test.okay.cds",
    "test.okay.mrna",
    "test.pubids",
    "test.pubids.old",
    "test.pubids.realt.log",
]


@env.task
async def test_evigene_tr2aacds_transcriptome() -> None:
    fasta = await fixture("genomics/homo_sapiens/genome/transcriptome.fasta")
    result = await evigene_tr2aacds(fasta, prefix="test")
    dropset = {f.name: f for f in await result.dropset.list_files()}
    assert sorted(dropset) == sorted(DROPSET), f"dropset holds {sorted(dropset)}"
    for name, md5 in DROPSET.items():
        await assert_md5(dropset[name], md5, label=name)
    okayset = sorted(f.name for f in await result.okayset.list_files())
    assert okayset == OKAYSET, f"okayset holds {okayset}"


tests = [test_evigene_tr2aacds_transcriptome]
