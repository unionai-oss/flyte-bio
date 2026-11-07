"""Tests for flyte_bio.modules.cellpose.

The case is upstream's cellpose test (modules commit 9a1cc9f, Cellpose 4.0.9),
with its ``ext.args`` ``--no_npy --save_flows``. As upstream, it checks the
mask files' names and that each is nonempty: masks from a neural network
aren't reproducible byte for byte. It runs on CPU; upstream's GPU test is the
same case with ``--use_gpu``.
"""

from flyte_bio.modules.cellpose import cellpose
from tests.framework import assert_nonempty, env, fixture

IMAGE = "imaging/segmentation/cycif_tonsil_registered.ome.tif"
MASKS = [
    "cycif_tonsil_registered.ome_cp_masks.tif",
    "cycif_tonsil_registered.ome_dP_cp_masks.tif",
    "cycif_tonsil_registered.ome_flows_cp_masks.tif",
]


@env.task
async def test_cellpose_with_flows_no_model() -> None:
    result = await cellpose(await fixture(IMAGE), args="--no_npy --save_flows")
    names = [f.name for f in result.masks]
    assert names == MASKS, f"got masks {names}, expected {MASKS}"
    for mask in result.masks:
        await assert_nonempty(mask, label=mask.name)
    assert result.cells == [], "wrote a _seg.npy despite --no_npy"


tests = [test_cellpose_with_flows_no_model]
