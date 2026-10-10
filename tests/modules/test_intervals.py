"""Tests for flyte_bio.modules.intervals.

Upstream has no module snapshots for these (its interval prep is local
awk), so the tests check the intervals themselves.
"""

from flyte.io import File

from flyte_bio.modules.intervals import build_intervals, create_intervals_bed
from tests.framework import env, fixture


async def text(f: File) -> str:
    async with f.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_build_intervals() -> None:
    fai = await fixture("genomics/sarscov2/genome/genome.fasta.fai")
    bed = await build_intervals(fai, prefix="genome")
    assert bed.path.endswith("/genome.bed"), bed.path
    assert await text(bed) == "MT192765.1\t0\t29829\n", await text(bed)


@env.task
async def test_create_intervals_bed() -> None:
    # Upstream's default test runs genome.multi_intervals.bed at --nucleotides_per_second 20:
    # chr22:1-15000 (~750 s) and chr22:20000-40001 (~1000 s) become two chunks, longest first.
    bed = await fixture("genomics/homo_sapiens/genome/genome.multi_intervals.bed")
    split = await create_intervals_bed(bed, nucleotides_per_second=20)
    names = [f.path.rsplit("/", 1)[-1] for f in split]
    assert names == ["chr22_20001-40001.bed", "chr22_2-15000.bed"], names
    assert [await text(f) for f in split] == ["chr22\t20000\t40001\n", "chr22\t1\t15000\n"]
    # At the default rate both are short, so they share one chunk.
    (single,) = await create_intervals_bed(bed)
    assert single.path.endswith("/chr22_2-15000.bed"), single.path
    assert await text(single) == "chr22\t1\t15000\nchr22\t20000\t40001\n", await text(single)


@env.task
async def test_create_intervals_from_interval_list() -> None:
    intervals = await fixture("genomics/homo_sapiens/genome/genome.interval_list")
    split = await create_intervals_bed(intervals)
    assert split, "no intervals from an interval_list"
    for f in split:
        for line in (await text(f)).splitlines():
            start, end = line.split("\t")[1:3]
            assert int(start) < int(end), line


tests = [test_build_intervals, test_create_intervals_bed, test_create_intervals_from_interval_list]
