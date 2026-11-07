"""Tests for flyte_bio.modules.hhsuite.

The cases and expected md5s are upstream's hhsuite/reformat tests (modules
commit 9a1cc9f, HH-suite 3.3.0), on the multiplesequencealign pipeline's
sea-toxin reference alignment; nf-test hashes ``.gz`` outputs decompressed.
The second case converts the first case's (gzipped) Stockholm output, as
upstream's setup step does.
"""

from flyte_bio.modules.hhsuite import hhsuite_reformat
from tests.framework import assert_gunzipped_md5, env, fixture

SETOXIN = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin.ref"


@env.task
async def test_hhsuite_reformat_setoxin_fas_sto() -> None:
    out = await hhsuite_reformat(await fixture(SETOXIN), "fas", "sto", prefix="test")
    assert out.name == "test.sto.gz", out.name
    await assert_gunzipped_md5(out, "41dcb53288298daebbdd912e74f0d018", label="test.sto.gz")


@env.task
async def test_hhsuite_reformat_setoxin_sto_a3m() -> None:
    sto = await hhsuite_reformat(await fixture(SETOXIN), "fas", "sto", prefix="test")
    out = await hhsuite_reformat(sto, "sto", "a3m", prefix="test")
    assert out.name == "test.a3m.gz", out.name
    await assert_gunzipped_md5(out, "c477b891bbb86399c62ede8ab8981aa5", label="test.a3m.gz")


tests = [test_hhsuite_reformat_setoxin_fas_sto, test_hhsuite_reformat_setoxin_sto_a3m]
