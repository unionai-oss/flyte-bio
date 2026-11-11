"""Tests for flyte_bio.modules.centrifuge.

The case and expected md5s are upstream's centrifuge/build test (modules
commit 9a1cc9f, Centrifuge 1.0.4.2): the SARS-CoV-2 genome indexed over the
test taxonomy.

The classification cases are upstream's centrifuge/centrifuge and
centrifuge/kreport tests (same commit), against upstream's test index
``delete_me/minigut_cf.tar.gz`` (extracted with :func:`fixture_dir`, as
upstream untars it), saving the aligned and unaligned reads. As upstream,
they check the outputs' names only.
"""

from flyte.io import Dir

from flyte_bio.modules.centrifuge import centrifuge_build, centrifuge_centrifuge, centrifuge_kreport
from tests.framework import assert_md5, env, fixture, fixture_dir

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


FASTQ = "genomics/sarscov2/illumina/fastq/"


async def minigut_index() -> Dir:
    return await fixture_dir("delete_me/minigut_cf.tar.gz", "minigut_cf")


@env.task
async def test_centrifuge_centrifuge_single_end() -> None:
    # upstream case: centrifuge/centrifuge "sarscov2_fastq_se"
    result = await centrifuge_centrifuge(
        await fixture(FASTQ + "test_1.fastq.gz"),
        await minigut_index(),
        save_unaligned=True,
        save_aligned=True,
        prefix="test",
    )
    assert (result.report.name, result.results.name) == ("test.report.txt", "test.results.txt")
    assert [f.name for f in result.mapped] == ["test.mapped.fastq.gz"], [f.name for f in result.mapped]
    assert [f.name for f in result.unmapped] == ["test.unmapped.fastq.gz"], [f.name for f in result.unmapped]


@env.task
async def test_centrifuge_centrifuge_paired_end() -> None:
    # upstream case: centrifuge/centrifuge "sarscov2_fastq_pe"
    result = await centrifuge_centrifuge(
        await fixture(FASTQ + "test_1.fastq.gz"),
        await minigut_index(),
        reads_2=await fixture(FASTQ + "test_2.fastq.gz"),
        save_unaligned=True,
        save_aligned=True,
        prefix="test",
    )
    assert (result.report.name, result.results.name) == ("test.report.txt", "test.results.txt")
    mapped = [f.name for f in result.mapped]
    unmapped = [f.name for f in result.unmapped]
    assert mapped == ["test.mapped.fastq.1.gz", "test.mapped.fastq.2.gz"], mapped
    assert unmapped == ["test.unmapped.fastq.1.gz", "test.unmapped.fastq.2.gz"], unmapped


@env.task
async def test_centrifuge_kreport_single_end() -> None:
    # upstream case: centrifuge/kreport "sarscov2_fastq_se"
    index = await minigut_index()
    result = await centrifuge_centrifuge(
        await fixture(FASTQ + "test_1.fastq.gz"), index, save_unaligned=True, save_aligned=True, prefix="test"
    )
    report = await centrifuge_kreport(result.results, index, prefix="test")
    assert report.name == "test.txt", report.name


tests = [
    test_centrifuge_build,
    test_centrifuge_centrifuge_single_end,
    test_centrifuge_centrifuge_paired_end,
    test_centrifuge_kreport_single_end,
]
