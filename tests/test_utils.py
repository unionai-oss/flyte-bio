"""Tests for the pipeline utilities: flyte_bio.samplesheet and flyte_bio.publish."""

import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.publish import publish
from flyte_bio.samplesheet import read_samplesheet
from tests.framework import env


async def write(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


async def tree(d: Dir) -> dict[str, str]:
    """Relative path -> content for every file under ``d``."""
    base = d.path.rstrip("/")
    out = {}
    async for f in d.walk():
        async with f.open("rb") as fh:
            out[f.path[len(base) :].lstrip("/")] = bytes(await fh.read()).decode()
    return out


@env.task
async def test_read_samplesheet() -> None:
    sheet = await write("samples.csv", "sample,fastq_1,fastq_2\n A ,s3://b/a_1.fq.gz, \nB,s3://b/b_1.fq.gz,s3://b/b_2.fq.gz\n")
    rows = await read_samplesheet(sheet)
    assert [dict(r) for r in rows] == [
        {"sample": "A", "fastq_1": "s3://b/a_1.fq.gz", "fastq_2": ""},
        {"sample": "B", "fastq_1": "s3://b/b_1.fq.gz", "fastq_2": "s3://b/b_2.fq.gz"},
    ], rows
    assert [r.line for r in rows] == [2, 3]
    try:
        rows[0]["strandedness"]
    except KeyError as e:
        assert "no column 'strandedness'" in str(e) and "line 2" in str(e), e
    else:
        raise AssertionError("a missing column should raise KeyError")


@env.task
async def test_publish() -> None:
    a = await write("a.txt", "a")
    root = Path(tempfile.mkdtemp())
    (root / "x").mkdir()
    (root / "x" / "1.txt").write_text("1")
    (root / "2.txt").write_text("2")
    d = await Dir.from_local(str(root))

    published = await publish({"top/renamed.txt": a, "dirs/d": d, "skipped": None})
    assert await tree(published) == {"top/renamed.txt": "a", "dirs/d/x/1.txt": "1", "dirs/d/2.txt": "2"}

    try:
        await publish({"same/2.txt": a, "same": d})
    except ValueError as e:
        assert "same/2.txt" in str(e), e
    else:
        raise AssertionError("two outputs at one path should raise ValueError")


tests = [test_read_samplesheet, test_publish]
