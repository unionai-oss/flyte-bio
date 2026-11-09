"""Tests for flyte_bio.modules.memote.

The cases are upstream's memote/run and memote/report tests (modules commit
9a1cc9f, MEMOTE 0.17.0), on upstream's gapseq-built model fixture. As
upstream, the outputs are checked for presence only (they record run times
and dates); the results are also checked to be gzipped JSON.
"""

import gzip
import json
import os

from flyte_bio.modules.memote import memote_report, memote_run
from tests.framework import assert_nonempty, env, fixture

MODEL = "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/gapseq/test.xml"


@env.task
async def test_memote_run() -> None:
    # upstream case: memote/run "candidatus_portiera - xml - json"
    result = await memote_run(await fixture(MODEL), prefix="test")
    assert os.path.basename(result.path) == "test.json.gz", f"unexpected name {result.path}"
    async with result.open("rb") as fh:
        report = json.loads(gzip.decompress(bytes(await fh.read())))
    assert "tests" in report, f"memote results lack 'tests': {sorted(report)}"


@env.task
async def test_memote_report() -> None:
    # upstream case: memote/report "candidatus_portiera - xml - html"
    result = await memote_report(await fixture(MODEL), prefix="test")
    assert os.path.basename(result.path) == "test.html", f"unexpected name {result.path}"
    await assert_nonempty(result, label="test.html")


tests = [test_memote_run, test_memote_report]
