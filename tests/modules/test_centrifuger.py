"""Tests for flyte_bio.modules.centrifuger.

The case is upstream's centrifuger/build test (modules commit 9a1cc9f,
Centrifuger 1.1.2): the SARS-CoV-2 genome indexed over the test taxonomy.
Upstream snapshots only the index's file names (its contents aren't
stable), so the test checks those.
"""

from flyte_bio.modules.centrifuger import centrifuger_build
from tests.framework import assert_nonempty, env, fixture

METAGENOME = "genomics/sarscov2/metagenome/"


@env.task
async def test_centrifuger_build() -> None:
    # upstream case: centrifuger/build "sarscov2 - fasta with reference file"
    index = await centrifuger_build(
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        taxonomy_nodes=await fixture(METAGENOME + "nodes.dmp"),
        taxonomy_names=await fixture(METAGENOME + "names.dmp"),
        conversion_table=await fixture(METAGENOME + "seqid2taxid.map"),
        prefix="test",
    )
    files = await index.list_files()
    names = sorted(f.name for f in files)
    assert names == ["test.1.cfr", "test.2.cfr", "test.3.cfr", "test.4.cfr"], f"centrifuger index: {names}"
    for file in files:
        await assert_nonempty(file, label=file.name)


tests = [test_centrifuger_build]
