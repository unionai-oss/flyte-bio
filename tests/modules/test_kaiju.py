"""Tests for flyte_bio.modules.kaiju.

The cases and expected md5s are upstream's kaiju/mkfmi tests (modules commit
9a1cc9f, Kaiju 1.10.0): the SARS-CoV-2 proteome indexed, without and with
the intermediate files. The test directory's config sets a protein alphabet
(``-a``), but the tests don't load it, and their md5s are those of the
default alphabet; so these tests pass no options either.
"""

from flyte.io import Dir

from flyte_bio.modules.kaiju import kaiju_mkfmi
from tests.framework import assert_md5, env, fixture

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


tests = [test_kaiju_mkfmi, test_kaiju_mkfmi_keep_intermediate]
