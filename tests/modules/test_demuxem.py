"""Tests for flyte_bio.modules.demuxem.

The case is upstream's demuxem test (modules commit 9a1cc9f), with its test
config's ``ext.args``. As upstream, only the outputs' names are checked
(demuxEM's zarr archives aren't byte-stable); here they must also be
nonempty.
"""

from flyte_bio.modules.demuxem import demuxem
from tests.framework import assert_nonempty, env, fixture

DATA = "genomics/homo_sapiens/10xgenomics/cellranger/hashing_demultiplexing/"
ARGS = (
    "--threads 1 --alpha-on-samples 0.0 --min-num-genes 100  --min-num-umis 100 --min-signal-hashtag 10.0"
    " --random-state 0"
)


@env.task
async def test_demuxem_h5_csv() -> None:
    # upstream case: demuxem "Standard_Multiome - h5 - csv"
    h5 = await fixture(DATA + "438-21-raw_feature_bc_matrix.h5")
    csv = await fixture(DATA + "438_21_raw_HTO.csv")
    result = await demuxem(h5, csv, prefix="sample1", args=ARGS)
    assert result.zarr.name == "sample1_demux.zarr.zip", result.zarr.name
    assert result.out_zarr.name == "sample1.out.demuxEM.zarr.zip", result.out_zarr.name
    await assert_nonempty(result.zarr, label=result.zarr.name)
    await assert_nonempty(result.out_zarr, label=result.out_zarr.name)


tests = [test_demuxem_h5_csv]
