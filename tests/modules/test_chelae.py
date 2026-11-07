"""Tests for flyte_bio.modules.chelae.

The cases and expected md5s are upstream's chelae/trim tests (modules commit
9a1cc9f, chelae 0.1.0). nf-test computes the reads' md5s decompressed; as
upstream, the JSON report (which varies between runs) is checked by name
only, and the reads of the plain- and BGZF-input cases by their record count.
The module's 4 CPUs match upstream's test CPUs.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.chelae import ChelaeTrimResult, chelae_trim, chelae_trim_cmd
from flyte_bio.modules.htslib import htslib_bgziptabix, tabix_bgzip
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
SINGLE_TSV = "9771a90cf76999433b4ce645fafc2151"
PAIRED_TSV = "d3db5b48d69b0680e2b08f199765679b"
PAIRED_READS = ["e3da014fbb9b428e952c62e8f0fb6402", "2ebae722295ea66d84075a3b042e2b42"]


async def check(result: ChelaeTrimResult, tsv: str, reads: list[str] | None, names: list[str]) -> None:
    assert [f.name for f in result.reads] == names, [f.name for f in result.reads]
    assert result.json.name == "test.chelae.json", result.json.name
    await assert_md5(result.metrics, tsv, label="test.chelae.tsv")
    for file, md5 in zip(result.reads, reads or [], strict=False):
        await assert_gunzipped_md5(file, md5, label=file.name)


async def records(file: File) -> int:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).count(b"\n") // 4


async def paired() -> tuple[File, File]:
    return await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")


@env.task
async def test_chelae_single_end() -> None:
    # upstream case: "sarscov2 - single_end - fastq"
    result = await chelae_trim(await fixture(FASTQ + "test_1.fastq.gz"), None, prefix="test")
    await check(result, SINGLE_TSV, ["4161df271f9bfcd25d5845a1e220dbec"], ["test.chelae.fastq.gz"])


@env.task
async def test_chelae_paired_end() -> None:
    # upstream case: "sarscov2 - paired_end - fastq"
    result = await chelae_trim(*await paired(), prefix="test")
    await check(result, PAIRED_TSV, PAIRED_READS, ["test_R1.chelae.fastq.gz", "test_R2.chelae.fastq.gz"])


@env.task
async def test_chelae_paired_end_adapter_fasta() -> None:
    # upstream case: "sarscov2 - paired_end - fastq - adapter_fasta"
    adapters = await fixture("delete_me/fastp/adapters.fasta")
    result = await chelae_trim(*await paired(), prefix="test", adapter_fasta=adapters)
    await check(result, PAIRED_TSV, PAIRED_READS, ["test_R1.chelae.fastq.gz", "test_R2.chelae.fastq.gz"])


@env.task
async def test_chelae_plain_fastq() -> None:
    # upstream case: "sarscov2 - single_end - plain fastq" (input decompressed with bgzip -d)
    plain, _ = await tabix_bgzip(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    assert plain.name == "test.fastq", plain.name
    result = await chelae_trim(plain, None, prefix="test")
    await check(result, SINGLE_TSV, None, ["test.chelae.fastq.gz"])
    assert await records(result.reads[0]) == 100


@env.task
async def test_chelae_bgzf_fastq() -> None:
    # upstream case: "sarscov2 - single_end - bgzf fastq" (input recompressed with bgzip)
    bgzf, _ = await htslib_bgziptabix(
        await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", out_ext="fastq", make_index=False
    )
    result = await chelae_trim(bgzf, None, prefix="test")
    await check(result, SINGLE_TSV, None, ["test.chelae.fastq.gz"])
    assert await records(result.reads[0]) == 100


@env.task
async def test_chelae_mismatched_inputs_fail() -> None:
    # upstream case: "sarscov2 - mismatched inputs and outputs - fails" (two reads, single-end outputs)
    r1, r2 = await paired()
    try:
        await chelae_trim_cmd(reads_1=r1, reads_2=r2, adapter_fasta=None, single_end=True, prefix="test", args="")
    except Exception as e:  # chelae rejects two inputs for one output
        print(f"chelae failed as expected: {e}")
    else:
        raise AssertionError("chelae trim succeeded with two inputs and one output")


tests = [
    test_chelae_single_end,
    test_chelae_paired_end,
    test_chelae_paired_end_adapter_fasta,
    test_chelae_plain_fastq,
    test_chelae_bgzf_fastq,
    test_chelae_mismatched_inputs_fail,
]
