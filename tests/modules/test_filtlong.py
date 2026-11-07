"""Tests for flyte_bio.modules.filtlong.

The cases and expected md5s are upstream's filtlong tests (modules commit
9a1cc9f, Filtlong 0.2.1), with upstream's test config: ``--min_length 10``
and the prefix ``test_lr``. nf-test computes the reads' md5 decompressed; the
logs' md5s are of Filtlong's report as written.
"""

from flyte.io import File

from flyte_bio.modules.filtlong import FiltlongResult, filtlong
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

NANOPORE = "genomics/sarscov2/nanopore/fastq/test.fastq.gz"
FASTQ = "genomics/sarscov2/illumina/fastq/"
ARGS = "--min_length 10"
READS_MD5 = "7567d853ada6ac142332619d0b541d76"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: FiltlongResult, log_md5: str) -> None:
    assert "Scoring long reads" in (await text(result.log)).splitlines(), "filtlong log lacks 'Scoring long reads'"
    await assert_gunzipped_md5(result.reads, READS_MD5, label="reads")
    await assert_md5(result.log, log_md5, label="log")


@env.task
async def test_filtlong_nanopore() -> None:
    result = await filtlong(await fixture(NANOPORE), prefix="test_lr", args=ARGS)
    await check(result, "99b1a8e5c6debbb0754872498f0614d7")


@env.task
async def test_filtlong_nanopore_illumina_single_end() -> None:
    short = [await fixture(FASTQ + "test_1.fastq.gz")]
    result = await filtlong(await fixture(NANOPORE), short, prefix="test_lr", args=ARGS)
    await check(result, "1540be73c90d4ca7bf23b26151cc84b7")


@env.task
async def test_filtlong_nanopore_illumina_paired_end() -> None:
    short = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    result = await filtlong(await fixture(NANOPORE), short, prefix="test_lr", args=ARGS)
    await check(result, "5dc53086954a5b01c38aa7945ea4a7af")


tests = [test_filtlong_nanopore, test_filtlong_nanopore_illumina_single_end, test_filtlong_nanopore_illumina_paired_end]
