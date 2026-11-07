"""Tests for flyte_bio.modules.ribodetector.

The case is upstream's ribodetector test (modules commit 9a1cc9f, ribodetector
0.3.3, CPU build), with its assertions: the log's sequence counts, each mate's
record count, and each mate's sorted read names. Upstream snapshots the 4159
names per mate; the test compares the md5 of the newline-joined sorted names
with that of upstream's lists.
"""

import gzip
import hashlib
import re

from flyte.io import File

from flyte_bio.modules.ribodetector import ribodetector
from tests.framework import env, fixture

FASTQ = "genomics/homo_sapiens/illumina/fastq/"
LOG_LINES = [
    "4159 sequences loaded",
    "Processed 4159 sequences in total",
    "Detected 4159 non-rRNA sequences",
    "Detected 0 rRNA sequences",
]
NAMES_MD5S = ["9b9304aac48e4e5a2bbeb5b3e0900732", "4f32485b4595d787a6980fc7c8e094fd"]


async def ribodetector_read_names(fastq: File) -> list[str]:
    async with fastq.open("rb") as fh:
        text = gzip.decompress(bytes(await fh.read())).decode()
    lines = text.splitlines()
    return [line[1:].split()[0] for line in lines[0::4]]


@env.task
async def test_ribodetector_rnaseq_paired_end() -> None:
    # upstream case: ribodetector "ribodetector - rnaseq PE input"
    result = await ribodetector(
        await fixture(FASTQ + "test_rnaseq_1.fastq.gz"),
        length=150,
        reads_2=await fixture(FASTQ + "test_rnaseq_2.fastq.gz"),
        prefix="test",
    )
    async with result.log.open("rb") as fh:
        log = re.sub(r"\x1b\[[0-9;]*m", "", bytes(await fh.read()).decode())
    for line in LOG_LINES:
        assert line in log, f"ribodetector log lacks {line!r}"
    assert [f.name for f in result.reads] == ["test.nonrna.1.fastq.gz", "test.nonrna.2.fastq.gz"]
    for mate, (fastq, md5) in enumerate(zip(result.reads, NAMES_MD5S, strict=True), start=1):
        names = await ribodetector_read_names(fastq)
        assert len(names) == 4159, f"mate {mate}: {len(names)} records, expected 4159"
        actual = hashlib.md5("\n".join(sorted(names)).encode()).hexdigest()
        assert actual == md5, f"mate {mate}: sorted read names md5 {actual}, expected {md5}"


tests = [test_ribodetector_rnaseq_paired_end]
