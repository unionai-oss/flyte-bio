"""Tests for flyte_bio.modules.trimgalore.

Expected md5s are the upstream trimgalore snapshot values, which nf-test
computes over the decompressed reads.
"""

from flyte_bio.modules.trimgalore import trimgalore
from tests.framework import assert_gunzipped_md5, assert_nonempty, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_trimgalore_single_end() -> None:
    # upstream case: trimgalore "sarscov2 - fastq - single-end"
    result = await trimgalore(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    await assert_gunzipped_md5(result.reads_1, "566d44cca0d22c522d6cf0e50c7165dc", label="single-end trimmed")
    assert result.reads_2 is None
    await assert_nonempty(result.reports[0], label="single-end report")
    assert await result.reads_after_filtering() > 0


@env.task
async def test_trimgalore_paired_end() -> None:
    # upstream case: trimgalore "sarscov2 - fastq - paired-end"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    result = await trimgalore(r1, r2, prefix="test")
    await assert_gunzipped_md5(result.reads_1, "566d44cca0d22c522d6cf0e50c7165dc", label="paired-end val_1")
    assert result.reads_2 is not None
    await assert_gunzipped_md5(result.reads_2, "3c023e8e890b897821df3dc98f48c2b3", label="paired-end val_2")
    assert len(result.reports) == 2
    assert await result.reads_after_filtering() > 0


@env.task
async def test_trimgalore_fastqc() -> None:
    # The pipeline runs FastQC on the trimmed reads via --fastqc_args.
    result = await trimgalore(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", fastqc=True)
    for name in ("test_trimmed_fastqc.html", "test_trimmed_fastqc.zip"):
        assert await result.results.get_file(name) is not None, f"trimgalore --fastqc_args wrote no {name}"


tests = [test_trimgalore_single_end, test_trimgalore_paired_end, test_trimgalore_fastqc]
