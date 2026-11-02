"""Tests for flyte_bio.modules.cellbender.

Upstream's test: a raw PBMC matrix trained for 5 epochs (``--epochs 5``),
on the CPU and, in a separate test, the GPU. CellBender's outputs aren't
reproducible, so upstream snapshots only their names; these tests check each
one exists and isn't empty.
"""

from flyte_bio.modules.cellbender import CellbenderResult, cellbender_removebackground
from tests.framework import assert_nonempty, env, fixture

RAW = "genomics/homo_sapiens/scrnaseq/h5ad/SRR28679759_raw_matrix.h5ad"


async def check(result: CellbenderResult, prefix: str) -> None:
    assert result.h5.path.endswith(f"/{prefix}.h5"), result.h5.path
    assert result.report is not None, "no report"
    for f in (
        result.h5,
        result.filtered_h5,
        result.posterior_h5,
        result.barcodes,
        result.metrics,
        result.report,
        result.pdf,
        result.log,
        result.checkpoint,
    ):
        await assert_nonempty(f, label=f.path)


@env.task
async def test_cellbender_removebackground() -> None:
    result = await cellbender_removebackground(await fixture(RAW), prefix="test", args="--epochs 5")
    await check(result, "test")


@env.task
async def test_cellbender_removebackground_gpu() -> None:
    result = await cellbender_removebackground(await fixture(RAW), prefix="test_gpu", gpu=True, args="--epochs 5")
    await check(result, "test_gpu")


tests = [test_cellbender_removebackground, test_cellbender_removebackground_gpu]
