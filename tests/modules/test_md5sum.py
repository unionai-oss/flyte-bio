"""Tests for flyte_bio.modules.md5sum.

The cases and expected md5s are upstream's md5sum tests (modules commit
9a1cc9f, coreutils 9.5). The BSD-style case passes upstream's ``--tag``.
"""

from flyte.io import File

from flyte_bio.modules.md5sum import file_name, md5sum
from tests.framework import assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


async def check_checksums(checksums: list[File], expected: dict[str, str]) -> None:
    names = [file_name(f) for f in checksums]
    assert names == list(expected), f"md5sum wrote {names}, expected {list(expected)}"
    for file, (name, md5) in zip(checksums, expected.items(), strict=True):
        await assert_md5(file, md5, label=name)


@env.task
async def test_md5sum_hello() -> None:
    result = await md5sum([await fixture("generic/txt/hello.txt")], as_separate_files=True, prefix="hello")
    await check_checksums(result.checksums, {"hello.txt.md5": "5c18e1db5460fb32fed66966483165fd"})


@env.task
async def test_md5sum_hello_bsd_style() -> None:
    result = await md5sum(
        [await fixture("generic/txt/hello.txt")], as_separate_files=True, prefix="hello", args="--tag"
    )
    await check_checksums(result.checksums, {"hello.txt.md5": "152a03f5dc7aa8db6612f63154ecbca2"})


@env.task
async def test_md5sum_paired_fastq_combined() -> None:
    files = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    result = await md5sum(files, as_separate_files=False, prefix="test")
    await check_checksums(result.checksums, {"test.md5": "dfb98e45cbb77a9a63ae7029aee38bd1"})


@env.task
async def test_md5sum_paired_fastq_separate() -> None:
    files = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    result = await md5sum(files, as_separate_files=True, prefix="test")
    await check_checksums(
        result.checksums,
        {
            "test_1.fastq.gz.md5": "63a3f187d06ae5cffe77fcd39bbdd63a",
            "test_2.fastq.gz.md5": "055dc46636836dbd1d9ede2e6ce8cd4e",
        },
    )


tests = [
    test_md5sum_hello,
    test_md5sum_hello_bsd_style,
    test_md5sum_paired_fastq_combined,
    test_md5sum_paired_fastq_separate,
]
