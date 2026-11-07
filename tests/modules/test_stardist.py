"""Tests for flyte_bio.modules.stardist.

The case and expected md5 are upstream's stardist test (modules commit
9a1cc9f, StarDist 0.9.2, TensorFlow 2.20.0): the pretrained
``2D_versatile_fluo`` model on a nuclear image. It runs on CPU, as upstream's
test does.
"""

from flyte_bio.modules.stardist import stardist
from tests.framework import assert_md5, env, fixture


@env.task
async def test_stardist_2d_tif() -> None:
    mask = await stardist(await fixture("imaging/segmentation/nuclear_image.tif"), model="2D_versatile_fluo")
    assert mask.name == "nuclear_image.stardist.tif", f"got {mask.name}"
    await assert_md5(mask, "ca29a05346d9e0d4c15bb133122560ce", label=mask.name)


tests = [test_stardist_2d_tif]
