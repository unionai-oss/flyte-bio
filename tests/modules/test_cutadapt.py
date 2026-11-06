"""Tests for flyte_bio.modules.cutadapt.

The cases and expected md5s are upstream's cutadapt tests (modules commit
9a1cc9f, cutadapt 5.2). nf-test computes the md5s over the decompressed reads.
Upstream's test ``nextflow.config`` sets ``-q 25``, but the tests don't load
it (its snapshots have no quality-trimmed line), so these run without args.
As upstream, the log is checked from its summary on, skipping the lines that
record the Python version and the core count.
"""

from flyte.io import File

from flyte_bio.modules.cutadapt import CutadaptResult, cutadapt
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
TRIMMED_1_MD5 = "4161df271f9bfcd25d5845a1e220dbec"
TRIMMED_2_MD5 = "2ebae722295ea66d84075a3b042e2b42"


async def check_log(result: CutadaptResult, expected: list[str]) -> None:
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    block = "\n".join(expected)
    assert block in log, f"cutadapt log lacks the summary {expected!r}:\n{log}"


@env.task
async def test_cutadapt_single_end() -> None:
    result = await cutadapt(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    assert len(result.reads) == 1, f"expected one trimmed file, got {len(result.reads)}"
    await assert_gunzipped_md5(result.reads[0], TRIMMED_1_MD5, label="test.trim.fastq.gz")
    await check_log(
        result,
        [
            "",
            "=== Summary ===",
            "",
            "Total reads processed:                     100",
            "Reads written (passing filters):           100 (100.0%)",
            "",
            "Total basepairs processed:        13,897 bp",
            "Total written (filtered):         13,897 bp (100.0%)",
        ],
    )


@env.task
async def test_cutadapt_paired_end() -> None:
    reads: list[File] = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    result = await cutadapt(reads[0], reads[1], prefix="test")
    assert len(result.reads) == 2, f"expected two trimmed files, got {len(result.reads)}"
    await assert_gunzipped_md5(result.reads[0], TRIMMED_1_MD5, label="test_1.trim.fastq.gz")
    await assert_gunzipped_md5(result.reads[1], TRIMMED_2_MD5, label="test_2.trim.fastq.gz")
    await check_log(
        result,
        [
            "",
            "=== Summary ===",
            "",
            "Total read pairs processed:                100",
            "Pairs written (passing filters):           100 (100.0%)",
            "",
            "Total basepairs processed:        27,645 bp",
            "  Read 1:        13,897 bp",
            "  Read 2:        13,748 bp",
            "Total written (filtered):         27,645 bp (100.0%)",
            "  Read 1:        13,897 bp",
            "  Read 2:        13,748 bp",
        ],
    )


tests = [test_cutadapt_single_end, test_cutadapt_paired_end]
