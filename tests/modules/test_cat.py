"""Tests for flyte_bio.modules.cat.

Expected md5s are the upstream cat/fastq and cat/cat snapshot values; for
``.gz`` outputs nf-test computes them over the decompressed content. The
cases join files in a fixed order, so a reordering bug changes the md5.

The cat/cat cases are those of the module's last working version: upstream
has since commented them out (the module now only fails with a deprecation
notice), but its snapshot still holds their md5s.

The concatenate cases are upstream find/concatenate's. It joins in name
order, so its gff3 + maf md5 differs from cat/cat's, which keeps list order.
"""

from flyte_bio.modules.cat import cat_cat, cat_fastq, concatenate
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
GENOME = "genomics/sarscov2/genome/"


@env.task
async def test_cat_fastq_single_end() -> None:
    # upstream case: cat/fastq "test_cat_fastq_single_end"
    shards = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    merged = await cat_fastq(shards)
    await assert_gunzipped_md5(merged, "ee314a9bd568d06617171b0c85f508da", label="single_end")


@env.task
async def test_cat_fastq_paired_end() -> None:
    # upstream case: cat/fastq "test_cat_fastq_paired_end" — R1 and R2 shards merged separately
    r1 = await cat_fastq([await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test2_1.fastq.gz")])
    r2 = await cat_fastq([await fixture(FASTQ + "test_2.fastq.gz"), await fixture(FASTQ + "test2_2.fastq.gz")])
    await assert_gunzipped_md5(r1, "3ad9406595fafec8172368f9cd0b6a22", label="paired_end R1")
    await assert_gunzipped_md5(r2, "a52cab0b840c7178b0ea83df1fdbe8d5", label="paired_end R2")


async def fasta_and_sizes():
    return [await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.sizes")]


async def gff3_and_maf():
    # In this order, which isn't alphabetical, so the join must keep list order.
    return [await fixture(GENOME + "genome.gff3.gz"), await fixture(GENOME + "alignment/last/contigs.genome.maf.gz")]


@env.task
async def test_cat_cat_unzipped_to_unzipped() -> None:
    # upstream case: cat/cat "sarscov2 - [ fasta, sizes ] - unzipped"
    out = await cat_cat(await fasta_and_sizes(), sample="test")
    await assert_md5(out, "f44b33a0e441ad58b2d3700270e2dbe2", label="test.fasta")


@env.task
async def test_cat_cat_zipped_to_zipped() -> None:
    # upstream case: cat/cat "sarscov2 - [ gff3_gz, maf_gz ] - zipped"
    out = await cat_cat(await gff3_and_maf(), sample="test")
    await assert_gunzipped_md5(out, "c439d3b60e7bc03e8802a451a0d9a5d9", label="test.gff3.gz")


@env.task
async def test_cat_cat_zipped_to_unzipped() -> None:
    # upstream case: cat/cat "sarscov2 - [ gff3_gz, maf_gz ] - unzipped" (ext.prefix cat.txt)
    out = await cat_cat(await gff3_and_maf(), sample="test", filename="cat.txt")
    await assert_md5(out, "c439d3b60e7bc03e8802a451a0d9a5d9", label="cat.txt")


@env.task
async def test_cat_cat_unzipped_to_zipped() -> None:
    # upstream case: cat/cat "sarscov2 - [ fasta, sizes ] - zipped" (ext.prefix cat.txt.gz)
    out = await cat_cat(await fasta_and_sizes(), sample="test", filename="cat.txt.gz")
    await assert_gunzipped_md5(out, "f44b33a0e441ad58b2d3700270e2dbe2", label="cat.txt.gz")


@env.task
async def test_cat_cat_single_file_zipped() -> None:
    # upstream case: cat/cat "sarscov2 - fasta - zipped" (ext.prefix cat.txt.gz)
    out = await cat_cat([await fixture(GENOME + "genome.fasta")], sample="test", filename="cat.txt.gz")
    await assert_gunzipped_md5(out, "6e9fe4042a72f2345f644f239272b7e6", label="cat.txt.gz")


@env.task
async def test_cat_cat_name_conflict() -> None:
    # upstream case: cat/cat "sarscov2 - genome - error: name conflict" — the default
    # output name, genome.fasta, is an input's name
    try:
        await cat_cat(await fasta_and_sizes(), sample="genome")
    except ValueError as e:
        if "can't be the same as for the output prefix" not in str(e):
            raise AssertionError(f"name conflict: unexpected error {e}") from e
        return
    raise AssertionError("name conflict: cat_cat accepted an output named like an input")


@env.task
async def test_concatenate_mixed_input() -> None:
    # upstream case: find/concatenate "test_mixed_input" — a gzipped and a plain input
    files = [await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "alignment/last/contigs.genome.maf.gz")]
    try:
        await concatenate(files, sample="test")
    except ValueError as e:
        if "A mix of both is not allowed." not in str(e):
            raise AssertionError(f"mixed input: unexpected error {e}") from e
        return
    raise AssertionError("mixed input: concatenate accepted gzipped and plain inputs together")


@env.task
async def test_concatenate_unzipped_unzipped() -> None:
    # upstream case: find/concatenate "test_cat_unzipped_unzipped"
    out = await concatenate(await fasta_and_sizes(), sample="test")
    await assert_md5(out, "f44b33a0e441ad58b2d3700270e2dbe2", label="test.fasta")


@env.task
async def test_concatenate_zipped_zipped() -> None:
    # upstream case: find/concatenate "test_cat_zipped_zipped" — named after the first input,
    # joined in name order (the maf first)
    out = await concatenate(await gff3_and_maf(), sample="test")
    assert out.path.endswith("/test.gff3.gz"), out.path
    await assert_gunzipped_md5(out, "1a983f97a01107efe42c8c3e7728b4c1", label="test.gff3.gz")


@env.task
async def test_concatenate_zipped_unzipped() -> None:
    # upstream case: find/concatenate "test_cat_zipped_unzipped" (ext.prefix cat.txt)
    out = await concatenate(await gff3_and_maf(), sample="test", filename="cat.txt")
    await assert_md5(out, "1a983f97a01107efe42c8c3e7728b4c1", label="cat.txt")


@env.task
async def test_concatenate_unzipped_zipped() -> None:
    # upstream case: find/concatenate "test_cat_unzipped_zipped" (ext.prefix cat.txt.gz)
    out = await concatenate(await fasta_and_sizes(), sample="test", filename="cat.txt.gz")
    await assert_gunzipped_md5(out, "f44b33a0e441ad58b2d3700270e2dbe2", label="cat.txt.gz")


@env.task
async def test_concatenate_one_file_unzipped_zipped() -> None:
    # upstream case: find/concatenate "test_cat_one_file_unzipped_zipped" (ext.prefix cat.txt.gz)
    out = await concatenate([await fixture(GENOME + "genome.fasta")], sample="test", filename="cat.txt.gz")
    await assert_gunzipped_md5(out, "6e9fe4042a72f2345f644f239272b7e6", label="cat.txt.gz")


tests = [
    test_cat_fastq_single_end,
    test_cat_fastq_paired_end,
    test_cat_cat_unzipped_to_unzipped,
    test_cat_cat_zipped_to_zipped,
    test_cat_cat_zipped_to_unzipped,
    test_cat_cat_unzipped_to_zipped,
    test_cat_cat_single_file_zipped,
    test_cat_cat_name_conflict,
    test_concatenate_mixed_input,
    test_concatenate_unzipped_unzipped,
    test_concatenate_zipped_zipped,
    test_concatenate_zipped_unzipped,
    test_concatenate_unzipped_zipped,
    test_concatenate_one_file_unzipped_zipped,
]
