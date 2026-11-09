"""Tests for flyte_bio.modules.muscle5.

The cases and expected md5s are upstream's muscle5/super5 tests (modules
commit 9a1cc9f, MUSCLE 5.1) on the seatoxin sequences; nf-test hashes ``.gz``
outputs decompressed. They rely on the module's default 4 CPUs, as upstream's
tests run with. The last case sets ``-perm all``, as upstream's
``perm_all.config``.
"""

from flyte_bio.modules.muscle5 import muscle5_super5
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

SEQUENCES = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin-ref.fa"
ALIGNMENT_MD5 = "46ba556df08f7aabbe5e1ba31d226b6a"


@env.task
async def test_muscle5_super5_seatoxin_uncompressed() -> None:
    out = await muscle5_super5(await fixture(SEQUENCES), prefix="test")
    assert [f.name for f in out] == ["test.aln"], f"outputs {[f.name for f in out]}"
    await assert_md5(out[0], ALIGNMENT_MD5, label="test.aln")


@env.task
async def test_muscle5_super5_seatoxin_compressed() -> None:
    out = await muscle5_super5(await fixture(SEQUENCES), prefix="test", compress=True)
    assert [f.name for f in out] == ["test.aln.gz"], f"outputs {[f.name for f in out]}"
    await assert_gunzipped_md5(out[0], ALIGNMENT_MD5, label="test.aln.gz")


@env.task
async def test_muscle5_super5_seatoxin_compressed_perm_all() -> None:
    out = await muscle5_super5(await fixture(SEQUENCES), prefix="test", compress=True, args="-perm all")
    expected = ["testabc.0.aln.gz", "testacb.0.aln.gz", "testbca.0.aln.gz", "testnone.0.aln.gz"]
    assert [f.name for f in out] == expected, f"outputs {[f.name for f in out]}"
    for f in out:
        await assert_gunzipped_md5(f, ALIGNMENT_MD5, label=f.name)


tests = [
    test_muscle5_super5_seatoxin_uncompressed,
    test_muscle5_super5_seatoxin_compressed,
    test_muscle5_super5_seatoxin_compressed_perm_all,
]
