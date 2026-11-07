"""Tests for flyte_bio.modules.fastqe.

The cases and expected md5s are upstream's fastqe tests (modules commit
9a1cc9f, fastqe 0.5.2).
"""

from flyte_bio.modules.fastqe import fastqe
from tests.framework import assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_fastqe_single_end() -> None:
    # upstream case: fastqe "sarscov2 single-end [fastq]"
    tsv = await fastqe([await fixture(FASTQ + "test_1.fastq.gz")], prefix="test")
    assert tsv.name == "test.tsv", tsv.name
    await assert_md5(tsv, "16d66ef19b5f2694aaeae9c3727c3a61", label="test.tsv")


@env.task
async def test_fastqe_paired_end() -> None:
    # upstream case: fastqe "sarscov2 paired-end [fastq]"
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    await assert_md5(await fastqe(reads, prefix="test"), "aceceddba6e4b5a968d54a6db46238da", label="test.tsv")


tests = [test_fastqe_single_end, test_fastqe_paired_end]
