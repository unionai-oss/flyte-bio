"""Tests for flyte_bio.modules.foldmason.

The cases and expected md5s are upstream's foldmason/createdb,
foldmason/easymsa and foldmason/msa2lddtreport tests (modules commit 9a1cc9f),
on the five seatoxin AlphaFold2 structures. nf-test hashes ``.gz`` outputs
decompressed. As upstream, the guide-tree case builds its tree with
:func:`famsa_guidetree`, and the report only has to exist (upstream snapshots
just its versions).
"""

import hashlib

from flyte.io import File

from flyte_bio.modules.famsa import famsa_guidetree
from flyte_bio.modules.foldmason import foldmason_createdb, foldmason_easymsa, foldmason_msa2lddtreport
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, env, fixture, fixture_dir

MSA_DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/"
STRUCTURES = MSA_DATA + "af2_structures/seatoxin-ref.tar.gz"
SEQUENCES = MSA_DATA + "setoxin-ref.fa"
REFERENCE_MSA = MSA_DATA + "setoxin.ref"
CREATEDB = {
    "test": "941321067eae439b4a8ccf6425d84751",
    "test.dbtype": "f1d3ff8443297732862df21dc4e57262",
    "test.index": "8d3c65fcda8b216fff2f1eb2c4dcc015",
    "test.lookup": "83047fddf9b6fbaa99f00cdf5a7dd274",
    "test.source": "7968a6cb6577ead230c31910380af948",
    "test_ca": "2738c516c1a15238a5c62c4ba6dee0b9",
    "test_ca.dbtype": "3fd85f9ee7ca8882c8caa747d0eef0b3",
    "test_ca.index": "cfdf544c3aa6d7e2034e4a01dac1d0ba",
    "test_h": "ab9ce99a99fc6ba6a98c4460410b6a16",
    "test_h.dbtype": "740bab4f9ec8808aedb68d6b1281aeb2",
    "test_h.index": "dc7c33ddb6a3dc54ad033120ef4c9af4",
    "test_ss": "75d329b63c0383c3e43090ba89238e14",
    "test_ss.dbtype": "f1d3ff8443297732862df21dc4e57262",
    "test_ss.index": "8d3c65fcda8b216fff2f1eb2c4dcc015",
}
MSA_3DI_MD5 = "87b7cb32dc09bc44f5700485c4c04944"
MSA_AA_MD5 = "65cd0c4b10f37d6e17f7c52bd64a7bf5"


async def seatoxin_structures() -> list[File]:
    structures = await fixture_dir(STRUCTURES, "seatoxin-ref")
    return sorted([f async for f in structures.walk()], key=lambda f: f.name)


@env.task
async def test_foldmason_createdb_seatoxin() -> None:
    db = await foldmason_createdb(await seatoxin_structures(), prefix="test")
    actual = {}
    async for file in db.walk():
        async with file.open("rb") as fh:
            actual[file.name] = hashlib.md5(bytes(await fh.read())).hexdigest()
    assert actual == CREATEDB, f"database md5s differ: {actual}"


@env.task
async def test_foldmason_easymsa_seatoxin_uncompressed() -> None:
    result = await foldmason_easymsa(await seatoxin_structures(), prefix="test")
    assert result.msa_3di.name == "test_3di.fa" and result.msa_aa.name == "test_aa.fa"
    await assert_md5(result.msa_3di, MSA_3DI_MD5, label="test_3di.fa")
    await assert_md5(result.msa_aa, MSA_AA_MD5, label="test_aa.fa")


@env.task
async def test_foldmason_easymsa_seatoxin_compressed() -> None:
    result = await foldmason_easymsa(await seatoxin_structures(), prefix="test", compress=True)
    assert result.msa_3di.name == "test_3di.fa.gz" and result.msa_aa.name == "test_aa.fa.gz"
    await assert_gunzipped_md5(result.msa_3di, MSA_3DI_MD5, label="test_3di.fa.gz")
    await assert_gunzipped_md5(result.msa_aa, MSA_AA_MD5, label="test_aa.fa.gz")


@env.task
async def test_foldmason_easymsa_seatoxin_guide_tree() -> None:
    tree = await famsa_guidetree(await fixture(SEQUENCES), prefix="test_tree")
    result = await foldmason_easymsa(await seatoxin_structures(), tree=tree, prefix="test")
    await assert_md5(result.msa_3di, "05af06303d43b84ba540e8a85957f789", label="test_3di.fa")
    await assert_md5(result.msa_aa, "ce75ecef9390218147cac9da63202ef7", label="test_aa.fa")


@env.task
async def test_foldmason_msa2lddtreport_seatoxin() -> None:
    structures = await seatoxin_structures()
    db = await foldmason_createdb(structures, prefix="test")
    tree = await famsa_guidetree(await fixture(SEQUENCES), prefix="test_tree")
    report = await foldmason_msa2lddtreport(
        await fixture(REFERENCE_MSA), db, db_name="test", structures=structures, tree=tree, prefix="test"
    )
    assert report.name == "test.html", f"output named {report.name}"
    await assert_nonempty(report, label="test.html")


tests = [
    test_foldmason_createdb_seatoxin,
    test_foldmason_easymsa_seatoxin_uncompressed,
    test_foldmason_easymsa_seatoxin_compressed,
    test_foldmason_easymsa_seatoxin_guide_tree,
    test_foldmason_msa2lddtreport_seatoxin,
]
