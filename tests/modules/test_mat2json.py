"""Tests for flyte_bio.modules.mat2json.

The case and expected md5 are upstream's mat2json test (modules commit
9a1cc9f, mat2json 1.0.0).
"""

import os

from flyte_bio.modules.mat2json import mat2json
from tests.framework import assert_md5, env, fixture

LSMQUANT = "https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/lsmquant/test_data/"


@env.task
async def test_mat2json_thresholds() -> None:
    converted = await mat2json(await fixture(LSMQUANT + "mat2json/thresholds.mat"))
    names = [f.name or os.path.basename(f.path) for f in converted]
    assert names == ["thresholds.json"], f"got {names}, expected ['thresholds.json']"
    await assert_md5(converted[0], "eb096d865e9f6a8ff6496b79eb448a92", label="thresholds.json")


tests = [test_mat2json_thresholds]
