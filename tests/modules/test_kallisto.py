"""Tests for flyte_bio.modules.kallisto.

The cases are upstream's kallisto/index and kallisto/quant tests (modules
commit 9a1cc9f, kallisto 0.52.0).

- index: upstream marks the index unstable and snapshots only its name, so the
  test checks the name and that the index isn't empty.
- quant: upstream marks the results directory and run info unstable and
  snapshots the results' file names and the log's md5; the strandedness case
  also checks that the run info's command line has ``--fr-stranded``.
"""

import json

from flyte.io import File

from flyte_bio.modules.kallisto import KallistoQuantResult, kallisto_index, kallisto_quant
from tests.framework import assert_md5, assert_nonempty, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
RESULT_FILES = ["abundance.h5", "abundance.tsv", "kallisto_quant.log", "run_info.json"]


@env.task
async def test_kallisto_index_transcriptome() -> None:
    # upstream case: kallisto/index "sarscov2 transcriptome.fasta"
    index = await kallisto_index(await fixture("genomics/sarscov2/genome/transcriptome.fasta"))
    assert index.name == "kallisto", index.name
    await assert_nonempty(index, label="kallisto index")


async def kallisto_quant_index() -> File:
    return await kallisto_index(await fixture("genomics/sarscov2/genome/transcriptome.fasta"))


async def check_kallisto_quant(result: KallistoQuantResult, log_md5: str) -> None:
    names = sorted(f.name for f in await result.results.list_files())
    assert names == RESULT_FILES, f"kallisto quant results hold {names}"
    assert result.json_info.name == "test.run_info.json", result.json_info.name
    assert result.log.name == "test.log", result.log.name
    await assert_md5(result.log, log_md5, label="log")


@env.task
async def test_kallisto_quant_single_end() -> None:
    # upstream case: kallisto/quant "sarscov2 single-end"
    result = await kallisto_quant(
        await fixture(FASTQ + "test_1.fastq.gz"),
        await kallisto_quant_index(),
        fragment_length=150,
        fragment_length_sd=75,
        prefix="test",
    )
    await check_kallisto_quant(result, "0cad96e4108b2ea1c152ad063198b21f")


@env.task
async def test_kallisto_quant_paired_end() -> None:
    # upstream case: kallisto/quant "sarscov2 paired-end"
    result = await kallisto_quant(
        await fixture(FASTQ + "test_1.fastq.gz"),
        await kallisto_quant_index(),
        reads_2=await fixture(FASTQ + "test_2.fastq.gz"),
        prefix="test",
    )
    await check_kallisto_quant(result, "cb46ef02463f63fec5aac31a4adcac94")


@env.task
async def test_kallisto_quant_single_end_stranded() -> None:
    # upstream case: kallisto/quant "sarscov2 single-end ext.args sets strandedness"
    result = await kallisto_quant(
        await fixture(FASTQ + "test_1.fastq.gz"),
        await kallisto_quant_index(),
        fragment_length=150,
        fragment_length_sd=75,
        prefix="test",
        args="--fr-stranded",
    )
    await check_kallisto_quant(result, "de11bd85ccf46173be2b144b8dc1f4f7")
    async with result.json_info.open("rb") as fh:
        call = json.loads(bytes(await fh.read()))["call"]
    assert "--fr-stranded" in call, f"run_info.json call lacks --fr-stranded: {call}"


tests = [
    test_kallisto_index_transcriptome,
    test_kallisto_quant_single_end,
    test_kallisto_quant_paired_end,
    test_kallisto_quant_single_end_stranded,
]
