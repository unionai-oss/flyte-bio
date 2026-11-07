"""Tests for flyte_bio.modules.mindagap.

The cases and expected md5s are upstream's mindagap/mindagap and
mindagap/duplicatefinder tests (modules commit 9a1cc9f, Mindagap 0.0.2), with
their ``ext.args``: ``3 40 --Xtilesize 2144`` and ``90``. The spot table comes
from upstream's molkart test data, as in its test.
"""

from flyte_bio.modules.mindagap import mindagap_duplicatefinder, mindagap_mindagap
from tests.framework import assert_md5, env, fixture

PANORAMA = "imaging/tiff/mindagap.mouse_heart.wga.tiff"
SPOTS = "https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt"


@env.task
async def test_mindagap_tiff() -> None:
    filled = await mindagap_mindagap(await fixture(PANORAMA), args="3 40 --Xtilesize 2144")
    assert filled.name == "mindagap.mouse_heart.wga_gridfilled.tiff", f"got {filled.name}"
    await assert_md5(filled, "310cf0017baa54af32176b43a5b0adfd", label=filled.name)


@env.task
async def test_mindagap_duplicatefinder_spots() -> None:
    marked = await mindagap_duplicatefinder(await fixture(SPOTS), args="90")
    assert marked.name == "spots_markedDups.txt", f"got {marked.name}"
    await assert_md5(marked, "4562caad05850d7dd7b6e9235e068a8b", label=marked.name)


tests = [test_mindagap_tiff, test_mindagap_duplicatefinder_spots]
