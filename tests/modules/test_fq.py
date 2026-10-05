"""Tests for flyte_bio.modules.fq.

Expected md5s are the upstream fq/subsample snapshot values, which nf-test
computes over the decompressed reads; the lint check mirrors upstream fq/lint.
"""

from flyte_bio.modules.fq import fq_lint, fq_subsample
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_fq_subsample_probability() -> None:
    # upstream case: fq/subsample "test_fq_subsample_probability"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    out_1, out_2 = await fq_subsample(r1, r2, args="-p 0.1 -s 123", prefix="test")
    assert out_2 is not None
    await assert_gunzipped_md5(out_1, "19326ff922a16c0cb81191f2a0a5c5fc", label="probability R1")
    await assert_gunzipped_md5(out_2, "ce7ff46296d89b68521ad55a3588bcfe", label="probability R2")


@env.task
async def test_fq_subsample_record_count() -> None:
    # upstream case: fq/subsample "test_fq_subsample_record_count"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    out_1, out_2 = await fq_subsample(r1, r2, args="-n 10 -s 123", prefix="test")
    assert out_2 is not None
    await assert_gunzipped_md5(out_1, "394c7a233f1c1c1a167a34cf2895d26d", label="record count R1")
    await assert_gunzipped_md5(out_2, "32724cbdb5ab954a0a659ebcd56ca422", label="record count R2")


@env.task
async def test_fq_subsample_single() -> None:
    # upstream case: fq/subsample "test_fq_subsample_single"
    out_1, out_2 = await fq_subsample(
        await fixture(FASTQ + "test_1.fastq.gz"), None, args="--probability 0.1 -s 123", prefix="test"
    )
    assert out_2 is None
    await assert_gunzipped_md5(out_1, "19326ff922a16c0cb81191f2a0a5c5fc", label="single")


@env.task
async def test_fq_lint_success() -> None:
    # upstream case: fq/lint "test_fq_lint_success"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    lint = await fq_lint(reads_1=r1, reads_2=[r2])
    async with lint.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    for marker in ("fq-lint start", "read 100 records", "fq-lint end"):
        assert marker in text, f"fq lint log lacks {marker!r}"


tests = [
    test_fq_subsample_probability,
    test_fq_subsample_record_count,
    test_fq_subsample_single,
    test_fq_lint_success,
]
