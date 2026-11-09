"""Tests for flyte_bio.modules.gnusort.

The cases and expected md5s are upstream's gnu/sort tests (modules commit
9a1cc9f, coreutils 9.5). Upstream's test sets ``params.gnu_sort_args``
(``-k1,1 -k2,2n``, or ``-t ';' -g -k 1,1 -k 2,2`` for the CSV) but never loads
the ``nextflow.config`` that would pass it on as ``ext.args``, so its
snapshots are of a plain ``sort`` (whole-line, byte order). The cases run
with no ``args`` to match: with the CSV's intended options, ``3;1;2`` would
sort before ``3;10;3``.
"""

from flyte_bio.modules.gnusort import gnu_sort
from tests.framework import assert_md5, env, fixture

UNSORTED = "generic/unsorted_data/unsorted_text/"


async def check(name: str, suffix: str, prefix: str, args: str, md5: str) -> None:
    result = await gnu_sort(await fixture(UNSORTED + name), suffix, prefix=prefix, args=args)
    assert result.name == f"{prefix}.{suffix}", result.name
    await assert_md5(result, md5, label=f"{prefix}.{suffix}")


@env.task
async def test_gnusort_genome() -> None:
    # upstream case: gnu/sort "unsorted_genome_sort"
    await check("test.genome", "genome.sorted", "genome_test", "", "fd97f7efafdbbfa71d9b560f10b4b048")


@env.task
async def test_gnusort_intervals() -> None:
    # upstream case: gnu/sort "unsorted_intervals_sort"
    await check("test.bed", "bed.sorted", "test", "", "abbce903ef263d38b2f71856387799ab")


@env.task
async def test_gnusort_csv() -> None:
    # upstream case: gnu/sort "unsorted_csv_sort"
    await check("test.csv", "csv.sorted", "test", "", "797450743dae151378fa5b7a9b5f38cc")


tests = [test_gnusort_genome, test_gnusort_intervals, test_gnusort_csv]
