"""Tests for flyte_bio.modules.pycoqc.

The case is upstream's pycoqc test (modules commit 9a1cc9f, pycoQC 2.5.2)
with its config's ``--min_pass_qual 0``: the HTML's name, and the JSON's
"All Reads" section against upstream's snapshot, which the test reads from
upstream's repository at that commit. Floats are compared within a relative
tolerance of 1e-9, since their last digits can differ between CPUs.
"""

import json
import math

from flyte_bio.modules.pycoqc import pycoqc
from tests.framework import env, fixture

SNAPSHOT = (
    "https://raw.githubusercontent.com/nf-core/modules/9a1cc9f2daeff56a0c53ee44cc583b577250782b/"
    "modules/nf-core/pycoqc/tests/main.nf.test.snap"
)


def pycoqc_close(actual, expected, where: str = "All Reads") -> None:
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and sorted(actual) == sorted(expected), f"{where}: keys {actual}"
        for key in expected:
            pycoqc_close(actual[key], expected[key], f"{where}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), f"{where}: length"
        for i, (a, e) in enumerate(zip(actual, expected, strict=True)):
            pycoqc_close(a, e, f"{where}[{i}]")
    elif isinstance(expected, float) or isinstance(actual, float):
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-12), f"{where}: {actual} != {expected}"
    else:
        assert actual == expected, f"{where}: {actual!r} != {expected!r}"


@env.task
async def test_pycoqc_sequencing_summary() -> None:
    # upstream case: pycoqc "sarscov2 nanopore [sequencing_summary.txt]"
    summary = await fixture("genomics/sarscov2/nanopore/sequencing_summary/test.sequencing_summary.txt")
    result = await pycoqc(summary, prefix="test", args="--min_pass_qual 0")
    assert result.html.name == "test.html", result.html.name
    async with result.json.open("rb") as fh:
        report = json.loads(bytes(await fh.read()))
    snapshot = await fixture(SNAPSHOT)
    async with snapshot.open("rb") as fh:
        expected = json.loads(bytes(await fh.read()))["all_reads"]["content"][0]
    pycoqc_close(report["All Reads"], expected)


tests = [test_pycoqc_sequencing_summary]
