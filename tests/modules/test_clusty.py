"""Tests for flyte_bio.modules.clusty.

The cases and expected md5s are upstream's clusty tests (modules commit
9a1cc9f, Clusty 1.2.2): a generic ANI table clustered alone, then with the
list of its objects. Each case sets column and threshold options in
``params.clusty_args``, which the test directory's ``nextflow.config`` makes
the ``ext.args``; but the tests don't load that config (no ``config`` line),
so their snapshots are of Clusty's default column and threshold choices, and
these tests pass no options either.
"""

from flyte_bio.modules.clusty import clusty
from tests.framework import assert_md5, env, fixture

ANI = "generic/tsv/ani.tsv"


@env.task
async def test_clusty_only_distances() -> None:
    # upstream case: clusty "generic - only distances"
    assignments = await clusty(await fixture(ANI), prefix="test")
    assert assignments.name == "test.tsv", assignments.name
    await assert_md5(assignments, "9c4d0da996fe7c4579b473f4148623ec", label="test.tsv")


@env.task
async def test_clusty_distances_objects() -> None:
    # upstream case: clusty "generic - distances - objects"
    objects = await fixture("generic/txt/ani_ids.txt")
    assignments = await clusty(await fixture(ANI), objects=objects, prefix="test")
    assert assignments.name == "test.tsv", assignments.name
    await assert_md5(assignments, "cea71211069e4aeed9e021ebebc038e7", label="test.tsv")


tests = [test_clusty_only_distances, test_clusty_distances_objects]
