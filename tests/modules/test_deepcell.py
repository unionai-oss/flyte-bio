"""Tests for flyte_bio.modules.deepcell.

The case and expected md5 are upstream's deepcell/mesmer test (modules commit
9a1cc9f, deepcell-applications 0.4.1, DeepCell 0.12.2), with its ``ext.prefix``
``mask`` and ``ext.args``. It runs on CPU, as upstream's (CPU-only) image does.
"""

from flyte_bio.modules.deepcell import deepcell_mesmer
from tests.framework import assert_md5, env, fixture

ARGS = "--image-mpp=0.65 --compartment=whole-cell --nuclear-channel 0 --membrane-channel 1"


@env.task
async def test_deepcell_mesmer_tif() -> None:
    image = await fixture("imaging/segmentation/cycif_tonsil_registered.ome.tif")
    mask = await deepcell_mesmer(image, prefix="mask", args=ARGS)
    await assert_md5(mask, "1550535389bd24d4ea4a8288502b0afa", label="mask.tif")


tests = [test_deepcell_mesmer_tif]
