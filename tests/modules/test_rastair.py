"""Tests for flyte_bio.modules.rastair.

The cases and expected md5s are upstream's rastair/call, rastair/mbias,
rastair/mbiasparser and rastair/methylkit tests (modules commit 9a1cc9f,
rastair 0.8.2). As upstream, the parser's test parses ``rastair mbias``'s
output and the methylKit test converts ``rastair call``'s (with ``--nOT 1
--nOB 1``); the PDF plot is only checked to exist (upstream marks it
unstable). Gzipped outputs are compared decompressed, as nf-test hashes them.
"""

from flyte.io import File

from flyte_bio.modules.rastair import rastair_call, rastair_mbias, rastair_mbiasparser, rastair_methylkit
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/test.paired_end.methylated.sorted.bam"
GENOME = "genomics/sarscov2/genome/genome.fasta"


async def rastair_inputs() -> dict[str, File]:
    return {
        "bam": await fixture(BAM),
        "bai": await fixture(BAM + ".bai"),
        "fasta": await fixture(GENOME),
        "fai": await fixture(GENOME + ".fai"),
    }


async def rastair_test_calls() -> File:
    return await rastair_call(**await rastair_inputs(), prefix="test", trim_ot="1", trim_ob="1")


@env.task
async def test_rastair_call() -> None:
    # upstream case: "sarscov2 paired-end methylated [bam] / sarscov2 genome [fasta]"
    await assert_md5(await rastair_test_calls(), "7508cfdd0045c688f82e4cb4f1d37373", label="rastair_call.txt")


@env.task
async def test_rastair_mbias() -> None:
    # upstream case: "sarscov2 paired-end methylated [bam] / sarscov2 genome [fasta]"
    txt = await rastair_mbias(**await rastair_inputs(), prefix="test")
    await assert_md5(txt, "de29e280a36531145f2a0d11d4d6d53c", label=txt.name)


@env.task
async def test_rastair_mbiasparser() -> None:
    # upstream case: "sarscov2 paired-end methylated [bam] / sarscov2 genome [fasta]"
    mbias = await rastair_mbias(**await rastair_inputs(), prefix="test")
    result = await rastair_mbiasparser(mbias, prefix="test")
    await assert_md5(result.csv, "e46ba9ffc77f03688559a3007e005d91", label=result.csv.name)
    assert (result.trim_ot, result.trim_ob) == ("0,0,0,0", "0,0,0,0"), (result.trim_ot, result.trim_ob)
    assert result.pdf is not None, "no test.rastair_mbias_processed.pdf"


@env.task
async def test_rastair_methylkit() -> None:
    # upstream case: "sarscov2 paired-end methylated [bam] / sarscov2 genome [fasta]"
    methylkit = await rastair_methylkit(await rastair_test_calls(), prefix="test")
    await assert_gunzipped_md5(methylkit, "5e241a791148be0af9dc48f23226ea39", label=methylkit.name)


tests = [test_rastair_call, test_rastair_mbias, test_rastair_mbiasparser, test_rastair_methylkit]
