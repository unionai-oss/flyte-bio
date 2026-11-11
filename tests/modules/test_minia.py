"""Tests for flyte_bio.modules.minia.

The cases are upstream's minia tests (modules commit 9a1cc9f, minia 3.2.6):
the outputs' names, and that the contigs contain ``AAA`` and the unitigs
``CCC``.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.minia import minia
from tests.framework import env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


async def check_minia(reads: list[File]) -> None:
    result = await minia(reads, prefix="test")
    names = [result.contigs.name, result.unitigs.name, result.h5.name, result.log.name]
    assert names == ["test.contigs.fa.gz", "test.unitigs.fa.gz", "test.h5", "test-minia.log"], names
    for file, motif in ((result.contigs, "AAA"), (result.unitigs, "CCC")):
        async with file.open("rb") as fh:
            text = gzip.decompress(bytes(await fh.read())).decode()
        assert motif in text, f"{file.name} lacks {motif}"


@env.task
async def test_minia_single_end() -> None:
    # upstream case: minia "sarscov2 - fastq - se"
    await check_minia([await fixture(FASTQ + "test_interleaved.fastq.gz")])


@env.task
async def test_minia_paired_end() -> None:
    # upstream case: minia "sarscov2 - fastq - pe"
    await check_minia([await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")])


tests = [test_minia_single_end, test_minia_paired_end]
