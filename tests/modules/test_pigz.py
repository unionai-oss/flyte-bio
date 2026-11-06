"""Tests for flyte_bio.modules.pigz.

The cases and md5s are upstream's pigz/compress and pigz/uncompress tests
(modules commit 9a1cc9f, pigz 2.8). nf-test hashes ``.gz`` outputs after
decompressing them.
"""

from flyte_bio.modules.pigz import pigz_compress, pigz_uncompress
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture


@env.task
async def test_pigz_compress_fasta() -> None:
    # upstream case: pigz/compress "sarscov2 - genome - fasta"
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    gz = await pigz_compress(fasta)
    assert gz.path.endswith("/genome.fasta.gz"), gz.path
    await assert_gunzipped_md5(gz, "6e9fe4042a72f2345f644f239272b7e6", label="genome.fasta.gz")


@env.task
async def test_pigz_uncompress_fastq() -> None:
    # upstream case: pigz/uncompress "Should run without failures"
    fastq = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    out = await pigz_uncompress(fastq)
    assert out.path.endswith("/test_1.fastq"), out.path
    await assert_md5(out, "4161df271f9bfcd25d5845a1e220dbec", label="test_1.fastq")


tests = [test_pigz_compress_fasta, test_pigz_uncompress_fastq]
