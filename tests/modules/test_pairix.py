"""Tests for flyte_bio.modules.pairix.

The case and expected md5s are upstream's pairix test (modules commit
9a1cc9f, pairix 0.3.7): the pairs pass through unchanged (nf-test hashes them
decompressed) and the ``.px2`` index is hashed as is.
"""

from flyte_bio.modules.pairix import pairix
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

PAIRS = (
    "https://raw.githubusercontent.com/4dn-dcic/pairix/master/samples/4dn.bsorted.chr21_22_only.nontriangle.pairs.gz"
)


@env.task
async def test_pairix() -> None:
    result = await pairix(await fixture(PAIRS))
    name = result.index.path.rstrip("/").rsplit("/", 1)[-1]
    assert name == "4dn.bsorted.chr21_22_only.nontriangle.pairs.gz.px2", name
    await assert_gunzipped_md5(result.pairs, "ee90dc99987ece492d598ed939d54561", label="pairs")
    await assert_md5(result.index, "a6e41cc7cff16fd15b5ee505549ec99a", label="px2")


tests = [test_pairix]
