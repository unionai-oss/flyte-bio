"""Tests for flyte_bio.modules.numorph.

Upstream's tests are those of modules commit 9a1cc9f (NuMorph images 1.0.0,
3D U-Net 1.0.9); each module's image stack is upstream's zip, passed as the
archive.

- intensity: upstream's case. As upstream, it checks the output files' names
  and counts (the ``.mat`` files record paths and dates, so aren't
  reproducible byte for byte); here each is also checked to be nonempty.
- resample: upstream has only a stub test (it says the data it needs is too
  large). This case runs the stub's input for real: the cropped stitched
  stack, with upstream's parameter file narrowed to the one marker the stack
  holds (``topro``, channel ``C01``; with two markers NuMorph fails looking for
  the second). It checks the resampled NIfTI's name, its dimensions (5 x 5 x 3
  at 25 um) and its md5, from a run of the same image.
- 3dunet: upstream's case, which upstream runs on a GPU; here it runs on CPU.
  As upstream, it checks the two output files' names; here also that the
  counts table is nonempty.
"""

import os
import struct
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.numorph import numorph_3dunet, numorph_intensity, numorph_resample
from tests.framework import assert_md5, assert_nonempty, env, fixture

LSMQUANT = "https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/lsmquant/test_data/"
PARAMETERS = LSMQUANT + "parameterfiles/params_TEST1_lsmquant.csv"
RAW_TILES = "https://zenodo.org/records/14916478/files/ctip2_topro.zip"
STITCHED = LSMQUANT + "unet_test_images.zip"
MODEL = "https://zenodo.org/records/16893708/files/075_121_model.h5"

VARIABLES = ["adj_params.mat", "path_table.mat", "thresholds.mat"]
SAMPLES = [
    f"intensity_adjustment/{name}"
    for name in [
        "ctip2_1_1_11.tif",
        "ctip2_1_2_11.tif",
        "ctip2_2_1_11.tif",
        "ctip2_2_2_11.tif",
        "flatfield_1.png",
        "flatfield_2.png",
        "tile_adj_1.png",
        "tile_adj_2.png",
        "topro_1_1_11.tif",
        "topro_1_2_11.tif",
        "topro_2_1_11.tif",
        "topro_2_2_11.tif",
        "y_adj_1.png",
        "y_adj_2.png",
    ]
]


def relative(file: File, root: str) -> str:
    return file.path.split(f"/{root}/", 1)[-1]


async def read_bytes(file: File) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read())


async def single_marker_parameters() -> File:
    """Upstream's parameter file with only the ``topro`` marker (channel ``C01``)."""
    lines = (await read_bytes(await fixture(PARAMETERS))).decode().splitlines()
    replaced = {"markers": "markers,topro", "channel_num": "channel_num,C01"}
    lines = [replaced.get(line.split(",", 1)[0], line) for line in lines]
    path = Path(tempfile.mkdtemp()) / "params_TEST1_topro.csv"
    path.write_text("\n".join(lines) + "\n")
    return await File.from_local(str(path))


@env.task
async def test_numorph_intensity() -> None:
    result = await numorph_intensity(await fixture(RAW_TILES), await fixture(PARAMETERS), "TEST1")
    variables = [relative(f, "variables") for f in result.variables]
    samples = [relative(f, "samples") for f in result.samples]
    assert variables == VARIABLES, f"got variables {variables}, expected {VARIABLES}"
    assert samples == SAMPLES, f"got samples {samples}, expected {SAMPLES}"
    for file in [*result.variables, *result.samples, result.nm_variables]:
        await assert_nonempty(file, label=file.path)


@env.task
async def test_numorph_resample() -> None:
    resampled = await numorph_resample(await fixture(STITCHED), await single_marker_parameters(), "stitched")
    names = [f.name or os.path.basename(f.path) for f in resampled]
    assert names == ["stitched_C1_topro_25.nii"], f"got {names}, expected ['stitched_C1_topro_25.nii']"
    header = await read_bytes(resampled[0])
    dims = struct.unpack("<4h", header[40:48])
    assert dims == (3, 5, 5, 3), f"resampled NIfTI has dims {dims}, expected (3, 5, 5, 3)"
    await assert_md5(resampled[0], "4abe8ff18918b95cf0bccc6d4730b764", label=names[0])


@env.task
async def test_numorph_3dunet() -> None:
    result = await numorph_3dunet(await fixture(STITCHED), await fixture(MODEL), "TEST1")
    names = sorted(f.name or os.path.basename(f.path) for f in await result.results.list_files())
    assert names == ["TEST1.csv", "TEST1_counts.csv"], f"got {names}, expected ['TEST1.csv', 'TEST1_counts.csv']"
    await assert_nonempty(result.counts, label="TEST1_counts.csv")


tests = [test_numorph_intensity, test_numorph_resample, test_numorph_3dunet]
