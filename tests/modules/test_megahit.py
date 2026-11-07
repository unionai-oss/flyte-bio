"""Tests for flyte_bio.modules.megahit.

The cases and expectations are upstream's megahit tests (modules commit
9a1cc9f, MEGAHIT 1.2.9). Upstream snapshots no md5s: it checks that the final
contigs have ``>k`` headers, that each kind of intermediate contigs is written
(its ``.each`` checks pass whenever the list is non-empty), and that the log
ends with "ALL DONE. Time elapsed". The module's 4 CPUs match upstream's test
CPUs.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.megahit import MegahitResult, megahit
from tests.framework import env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


async def raw(file: File) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read())


async def check(result: MegahitResult) -> None:
    assert b">k" in gzip.decompress(await raw(result.contigs)), "megahit contigs lack '>k' headers"
    for label, files in [
        ("k*.contigs.fa.gz", result.k_contigs),
        ("k*.addi.fa.gz", result.addi_contigs),
        ("k*.local.fa.gz", result.local_contigs),
        ("k*.final.contigs.fa.gz", result.kfinal_contigs),
    ]:
        assert files, f"megahit wrote no intermediate_contigs/{label}"
    last = (await raw(result.log)).decode().splitlines()[-1]
    assert "ALL DONE. Time elapsed" in last, f"megahit log's last line is {last!r}"


@env.task
async def test_megahit_single_end() -> None:
    await check(await megahit(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test"))


@env.task
async def test_megahit_paired_end() -> None:
    reads_1 = await fixture(FASTQ + "test_1.fastq.gz")
    reads_2 = await fixture(FASTQ + "test_2.fastq.gz")
    await check(await megahit(reads_1, reads_2, prefix="test"))


@env.task
async def test_megahit_paired_end_coassembly() -> None:
    reads_1 = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test2_1.fastq.gz")]
    reads_2 = [await fixture(FASTQ + "test_2.fastq.gz"), await fixture(FASTQ + "test2_2.fastq.gz")]
    await check(await megahit(reads_1, reads_2, prefix="test"))


tests = [test_megahit_single_end, test_megahit_paired_end, test_megahit_paired_end_coassembly]
