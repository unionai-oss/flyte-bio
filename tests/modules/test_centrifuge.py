"""Tests for flyte_bio.modules.centrifuge.

The case and expected md5s are upstream's centrifuge/build test (modules
commit 9a1cc9f, Centrifuge 1.0.4.2): the SARS-CoV-2 genome indexed over the
test taxonomy.
"""

from flyte_bio.modules.centrifuge import centrifuge_build
from tests.framework import assert_md5, env, fixture

METAGENOME = "genomics/sarscov2/metagenome/"
INDEX_MD5S = {
    "test.1.cf": "07d811cd4e350d56267183d2ac7023a5",
    "test.2.cf": "23551529b34af3585d35f0e36461d9df",
    "test.3.cf": "445a974b5d64518944cb539df4149850",
    "test.4.cf": "c609dc3c247d5f3ac5426079de7fbe36",
}


@env.task
async def test_centrifuge_build() -> None:
    # upstream case: centrifuge/build "sarscov2 - fasta"
    index = await centrifuge_build(
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        conversion_table=await fixture(METAGENOME + "seqid2taxid.map"),
        taxonomy_tree=await fixture(METAGENOME + "nodes.dmp"),
        name_table=await fixture(METAGENOME + "names.dmp"),
        prefix="test",
    )
    names = sorted(f.name for f in await index.list_files())
    assert names == sorted(INDEX_MD5S), f"centrifuge index: {names}"
    for name, md5 in INDEX_MD5S.items():
        file = await index.get_file(name)
        assert file is not None
        await assert_md5(file, md5, label=name)


tests = [test_centrifuge_build]
