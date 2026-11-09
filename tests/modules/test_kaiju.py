"""Tests for flyte_bio.modules.kaiju.

The cases and expected md5s are upstream's kaiju/mkfmi tests (modules commit
9a1cc9f, Kaiju 1.10.0): the SARS-CoV-2 proteome indexed, without and with
the intermediate files. The test directory's config sets a protein alphabet
(``-a``), but the tests don't load it, and their md5s are those of the
default alphabet; so these tests pass no options either.

The kaiju/kaiju, kaiju/kaiju2table and kaiju/kaiju2krona cases (same commit)
classify test_1.fastq.gz against upstream's prebuilt test database
(``kaiju.tar.gz``: ``proteins.fmi`` with its ``nodes.dmp`` and
``names.dmp``), extracted with :func:`fixture_dir` where upstream untars it.
Upstream's kaiju tests also stage a second ``nodes.dmp``/``names.dmp`` beside
the archive, but the module's lookup picks the archive's (and kaiju2table and
kaiju2krona get the archive alone), so these tests pass the archive's. Kaiju's
per-read output isn't snapshotted upstream, which checks its name and one
classified read, as these tests do. Upstream's "paired-end" case passes the
same single-end reads as its single-end case, and so does its port here.
"""

from flyte.io import Dir, File

from flyte_bio.modules.kaiju import kaiju, kaiju2krona, kaiju2table, kaiju_mkfmi
from tests.framework import assert_md5, env, fixture, fixture_dir

KAIJU_READS = "genomics/sarscov2/illumina/fastq/test_1.fastq.gz"
KAIJU_CLASSIFIED_READ = "C\tERR5069949.2257580\t2697049"

DB_MD5S = {
    "names.dmp": "cac3fc251ff837dc31ef4b1832553d48",
    "nodes.dmp": "cf7b63f1194c1eb7baa5d3e1b745a69b",
    "test.fmi": "54fd89f5e4eab61af30175e8aa389598",
}
INTERMEDIATE_MD5S = {
    "test.bwt": "fe639339371f5b4352017522f54d0c99",
    "test.sa": "7d4a6481ad5ffdc2058bb35331a64d4d",
}


async def kaiju_test_db(keep_intermediate: bool) -> Dir:
    return await kaiju_mkfmi(
        await fixture("genomics/sarscov2/genome/proteome.fasta"),
        nodes_dmp=await fixture("delete_me/kaiju/nodes.dmp"),
        names_dmp=await fixture("delete_me/kaiju/names.dmp"),
        keep_intermediate=keep_intermediate,
        prefix="test",
    )


async def check_kaiju_db(db: Dir, expected: dict[str, str]) -> None:
    names = sorted(f.name for f in await db.list_files())
    assert names == sorted(expected), f"kaiju db: {names}"
    mismatches = []
    for name, md5 in expected.items():
        file = await db.get_file(name)
        assert file is not None
        try:
            await assert_md5(file, md5, label=name)
        except AssertionError as error:
            mismatches.append(str(error))
    assert not mismatches, "; ".join(mismatches)


@env.task
async def test_kaiju_mkfmi() -> None:
    # upstream case: kaiju/mkfmi "sarscov2 - proteome - fasta"
    await check_kaiju_db(await kaiju_test_db(False), DB_MD5S)


@env.task
async def test_kaiju_mkfmi_keep_intermediate() -> None:
    # upstream case: kaiju/mkfmi "sarscov2 - proteome - keep intermediates"
    await check_kaiju_db(await kaiju_test_db(True), DB_MD5S | INTERMEDIATE_MD5S)


async def kaiju_prebuilt_db() -> Dir:
    return await fixture_dir("genomics/sarscov2/genome/db/kaiju.tar.gz", "kaiju")


async def kaiju_classify_test_reads() -> File:
    return await kaiju(await fixture(KAIJU_READS), await kaiju_prebuilt_db(), prefix="test")


async def check_kaiju_results(results: File) -> None:
    assert results.name == "test.tsv", results.name
    async with results.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    assert KAIJU_CLASSIFIED_READ in text, "test.tsv lacks the classified read ERR5069949.2257580"


@env.task
async def test_kaiju_single_end() -> None:
    # upstream case: kaiju/kaiju "sarscov2 - fastq - single-end"
    await check_kaiju_results(await kaiju_classify_test_reads())


@env.task
async def test_kaiju_paired_end() -> None:
    # upstream case: kaiju/kaiju "sarscov2 - fastq - paired-end" (single-end reads, as upstream)
    await check_kaiju_results(await kaiju_classify_test_reads())


@env.task
async def test_kaiju2table_species() -> None:
    # upstream case: kaiju/kaiju2table "sarscov2 - fastq - single-end"
    summary = await kaiju2table(await kaiju_classify_test_reads(), await kaiju_prebuilt_db(), "species", prefix="test")
    assert summary.name == "test.txt", summary.name
    await assert_md5(summary, "0d9f8fd36fcf2888296ae12632c5f0a8", label="test.txt")


@env.task
async def test_kaiju2krona() -> None:
    # upstream case: kaiju/kaiju2krona "sarscov2 - fastq - single-end"
    txt = await kaiju2krona(await kaiju_classify_test_reads(), await kaiju_prebuilt_db(), prefix="test")
    assert txt.name == "test.txt", txt.name
    await assert_md5(txt, "68b2309d37767e444193fa6cea7c0494", label="test.txt")


tests = [
    test_kaiju_mkfmi,
    test_kaiju_mkfmi_keep_intermediate,
    test_kaiju_single_end,
    test_kaiju_paired_end,
    test_kaiju2table_species,
    test_kaiju2krona,
]
