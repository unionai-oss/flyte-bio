"""Tests for flyte_bio.modules.gfaffix.

The case and expected md5s are upstream's gfaffix test (modules commit
9a1cc9f, GFAffix 0.2.1).
"""

from flyte_bio.modules.gfaffix import gfaffix
from tests.framework import assert_md5, env, fixture


@env.task
async def test_gfaffix_assembly() -> None:
    # upstream case: gfaffix "test-gfaffix"
    result = await gfaffix(await fixture("genomics/sarscov2/illumina/gfa/assembly.gfa"), prefix="test")
    assert result.gfa.path.endswith("/test.gfaffix.gfa"), result.gfa.path
    assert result.affixes.path.endswith("/test.affixes.txt"), result.affixes.path
    await assert_md5(result.gfa, "ae5591e4dd169f66970f0d14dcdd83e1", label="test.gfaffix.gfa")
    await assert_md5(result.affixes, "c3540baa386462552d7ebe4e4b7eda6a", label="test.affixes.txt")


tests = [test_gfaffix_assembly]
