"""Tests for flyte_bio.modules.backsub.

The cases and expected md5s are upstream's backsub tests (modules commit
9a1cc9f, backsub 0.5.2). Upstream's tests set ``module_args`` but load no
config that passes it to the module, so backsub runs without ``args``, as
here. The failure case uses upstream's inline markers table, whose
``TestBackground`` isn't a marker.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.backsub import backsub
from tests.framework import assert_md5, env, fixture

DATA = "imaging/background_subtraction/"
INVALID_MARKERS = """channel,cycle,marker_name,exposure,background,remove
1,1,DNA 1,100,,
2,1,Na/K ATPase,100,,
3,1,CD3,100,,
4,1,CD45RO,100,,
5,2,DNA 2,100,,TRUE
6,2,Antigen Ki67,100,TestBackground,
7,2,Pan-cytokeratin,100,CD3,
8,2,Aortic smooth muscle actin,100,CD45RO,
"""


@env.task
async def test_backsub_run() -> None:
    # upstream case: "Should run without failures"
    image = await fixture(DATA + "cycif_tonsil_registered.ome.tif")
    result = await backsub(image, await fixture(DATA + "markers.csv"), prefix="test")
    await assert_md5(result.image, "09d6238472b196e2b1d6e6281d4fdf79", label="test.ome.tif")
    await assert_md5(result.markers, "825e5202259a930cfc9e79f5dc0c02b0", label="test.csv")


@env.task
async def test_backsub_missing_background_fails() -> None:
    # upstream case: "Should fail if background not in marker_name column"
    path = Path(tempfile.mkdtemp()) / "markers_invalid.csv"
    path.write_text(INVALID_MARKERS)
    markers = await File.from_local(str(path))
    image = await fixture(DATA + "cycif_tonsil_registered.ome.tif")
    try:
        await backsub(image, markers, prefix="test")
    except Exception as e:  # the task fails: backsub rejects the unknown background marker
        print(f"backsub failed as expected: {e}")
    else:
        raise AssertionError("backsub accepted a background marker missing from marker_name")


tests = [test_backsub_run, test_backsub_missing_background_fails]
