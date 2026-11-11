"""Tests for flyte_bio.modules.seqkit.

The cases and expected md5s are upstream's seqkit/stats, seqkit/replace,
seqkit/seq, seqkit/rmdup and seqkit/split2 tests (modules commit 9a1cc9f,
seqkit 2.13.0), with each case's ``ext.args`` (seqkit/seq's
``ext.args2 = '-n'`` is ``gzip_args``). seqkit/split2's uncompressed case
gunzips its input in the test, where upstream runs its gunzip module.
Upstream's replace patterns are shell-quoted (``-p 'A'``); here they are bare.
nf-test hashes ``.gz`` outputs decompressed.

The seqkit/fq2fa and seqkit/grep cases are upstream's at the same commit.
seqkit/grep's pattern-file cases write upstream's ``chr22.header`` (``chr22``
and a newline); its ``with_pattern`` case, ``ext.args = '-p chr22'``, passes
the pattern as the dedicated ``pattern`` input.

The seqkit/head and seqkit/sliding cases are upstream's at the same commit,
seqkit/sliding's with its test config's ``-s 2 -W 5``. seqkit/sliding writes
uncompressed output even for the gzipped FASTQ, as upstream.
"""

import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.seqkit import (
    seqkit_fq2fa,
    seqkit_grep,
    seqkit_head,
    seqkit_replace,
    seqkit_rmdup,
    seqkit_seq,
    seqkit_sliding,
    seqkit_sort,
    seqkit_split2,
    seqkit_stats,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
GENOME = "genomics/sarscov2/genome/"


def name_of(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


@env.task
async def test_seqkit_stats_single_end() -> None:
    stats = await seqkit_stats([await fixture(FASTQ + "test_1.fastq.gz")], prefix="test")
    await assert_md5(stats, "76d5ab444de880ef443f358841e44cf9", label="test.tsv")


@env.task
async def test_seqkit_stats_paired_end() -> None:
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    stats = await seqkit_stats(reads, prefix="test")
    await assert_md5(stats, "2d96d9f7d85456090e138b81bf97bcd5", label="test.tsv")


@env.task
async def test_seqkit_stats_nanopore() -> None:
    stats = await seqkit_stats([await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz")], prefix="test")
    await assert_md5(stats, "388047b0263b8d0a3041d6d80b316782", label="test.tsv")


@env.task
async def test_seqkit_stats_genome_fasta() -> None:
    stats = await seqkit_stats([await fixture(GENOME + "genome.fasta")], prefix="test")
    await assert_md5(stats, "6f6f8b42abd183c6ab90346a311ca522", label="test.tsv")


@env.task
async def test_seqkit_stats_transcriptome_fasta() -> None:
    stats = await seqkit_stats([await fixture(GENOME + "transcriptome.fasta")], prefix="test")
    await assert_md5(stats, "c79d860b2b0379810e4df3e8359fc8ed", label="test.tsv")


@env.task
async def test_seqkit_replace_fasta_gz() -> None:
    out = await seqkit_replace(await fixture(GENOME + "proteome.fasta.gz"), prefix="test", args="-s -p A -r N")
    assert name_of(out) == "test.fasta.gz", f"expected test.fasta.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, "c40eaff961f6f2a48bb7e8fd156ed5d7", label="test.fasta.gz")


@env.task
async def test_seqkit_replace_uncompressed() -> None:
    out = await seqkit_replace(await fixture(GENOME + "genome.fasta"), prefix="test", args="-s -p A -r N")
    assert name_of(out) == "test.fasta", f"expected test.fasta, got {name_of(out)}"
    await assert_md5(out, "b1518908253a4997fcad98270751112e", label="test.fasta")


@env.task
async def test_seqkit_replace_custom_extension() -> None:
    genome = await fixture(GENOME + "genome.fasta")
    out = await seqkit_replace(genome, prefix="test", out_ext="fasta", args="-s -p T -r N")
    assert name_of(out) == "test.fasta", f"expected test.fasta, got {name_of(out)}"
    await assert_md5(out, "05d3294a62c72f5489f067c1da3c2f6c", label="test.fasta")


@env.task
async def test_seqkit_seq_genome_fasta() -> None:
    out = await seqkit_seq(await fixture(GENOME + "genome.fasta"), prefix="test", gzip_args="-n")
    assert name_of(out) == "test.fasta", f"expected test.fasta, got {name_of(out)}"
    await assert_md5(out, "483f4a5dfe60171c86ee9b7e6dff908b", label="test.fasta")


@env.task
async def test_seqkit_seq_genome_fasta_gz() -> None:
    out = await seqkit_seq(await fixture(GENOME + "genome.fasta.gz"), prefix="test", gzip_args="-n")
    assert name_of(out) == "test.fasta.gz", f"expected test.fasta.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, "483f4a5dfe60171c86ee9b7e6dff908b", label="test.fasta.gz")


@env.task
async def test_seqkit_seq_fastq_gz() -> None:
    out = await seqkit_seq(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", gzip_args="-n")
    assert name_of(out) == "test.fastq.gz", f"expected test.fastq.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, "4161df271f9bfcd25d5845a1e220dbec", label="test.fastq.gz")


@env.task
async def test_seqkit_seq_name_conflict_fails() -> None:
    try:
        await seqkit_seq(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test_1", gzip_args="-n")
    except ValueError as e:
        assert "Input and output names are the same" in str(e), str(e)
    else:
        raise AssertionError("seqkit_seq accepted an output name equal to its input's")


NO_DUPLICATES_LOG = "a41135cfe024baaf42f135583fe73f0d"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check_rmdup(fastx: File, expected_name: str, md5: str, log_md5: str, removed: int) -> None:
    result = await seqkit_rmdup(fastx, prefix="test")
    assert name_of(result.fastx) == expected_name, f"expected {expected_name}, got {name_of(result.fastx)}"
    if expected_name.endswith(".gz"):
        await assert_gunzipped_md5(result.fastx, md5, label=expected_name)
    else:
        await assert_md5(result.fastx, md5, label=expected_name)
    assert name_of(result.log) == "test.log", f"expected test.log, got {name_of(result.log)}"
    log = await text(result.log)
    assert f"{removed} duplicated records removed" in log, log
    await assert_md5(result.log, log_md5, label="test.log")


@env.task
async def test_seqkit_rmdup_genome_fasta() -> None:
    fasta = await fixture(GENOME + "genome.fasta")
    await check_rmdup(fasta, "test.fasta", "483f4a5dfe60171c86ee9b7e6dff908b", NO_DUPLICATES_LOG, 0)


@env.task
async def test_seqkit_rmdup_repeated_fasta() -> None:
    path = Path(tempfile.mkdtemp()) / "repeated.fasta"
    path.write_text(">A\nAGCTAGCTAGCT\n>B\nAGCTAGCTAGCT\n>A\nAGCTAGCTAGCT")
    fasta = await File.from_local(str(path))
    await check_rmdup(fasta, "test.fasta", "7510a742291241e7d7556bf720caf65c", "314c0aaef0f832a217a3f6ce3f8bc117", 1)


@env.task
async def test_seqkit_rmdup_genome_fasta_gz() -> None:
    fasta = await fixture(GENOME + "genome.fasta.gz")
    await check_rmdup(fasta, "test.fasta.gz", "483f4a5dfe60171c86ee9b7e6dff908b", NO_DUPLICATES_LOG, 0)


@env.task
async def test_seqkit_rmdup_fastq_gz() -> None:
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    await check_rmdup(reads, "test.fastq.gz", "4161df271f9bfcd25d5845a1e220dbec", NO_DUPLICATES_LOG, 0)


@env.task
async def test_seqkit_rmdup_name_conflict_fails() -> None:
    try:
        await seqkit_rmdup(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test_1")
    except ValueError as e:
        assert "Input and output names are the same" in str(e), str(e)
    else:
        raise AssertionError("seqkit_rmdup accepted an output name equal to its input's")


async def check_split2(parts: list[File], expected: dict[str, str]) -> None:
    names = [name_of(f) for f in parts]
    assert names == list(expected), f"expected parts {list(expected)}, got {names}"
    for part, (name, md5) in zip(parts, expected.items()):
        if name.endswith(".gz"):
            await assert_gunzipped_md5(part, md5, label=name)
        else:
            await assert_md5(part, md5, label=name)


async def paired_reads() -> list[File]:
    return [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]


SPLIT_BY_LENGTH_1 = {
    "test_1.part_001.fastq.gz": "ffc87eb34d91d23ec9095bd8609a6a70",
    "test_1.part_002.fastq.gz": "8d71a0abe239e05e5c57c4d27c799a1d",
}
SPLIT_BY_PART_1 = {
    "test_1.part_001.fastq.gz": "b3c4d28be7ea9b960fbf1cb452adb53c",
    "test_1.part_002.fastq.gz": "c134c70c5b70c6b3c65979448b38917e",
    "test_1.part_003.fastq.gz": "822a40283637e2715e77d1ed0ed5bd52",
}
SPLIT_BY_SIZE_1 = {
    "test_1.part_001.fastq.gz": "ecc4b1841cd94704bba742ea4dcd48b0",
    "test_1.part_002.fastq.gz": "b3de467f2b6ab0d14e1f6ce14932a411",
}


@env.task
async def test_seqkit_split2_single_end_length() -> None:
    parts = await seqkit_split2([await fixture(FASTQ + "test_1.fastq.gz")], args="--by-length 8K")
    await check_split2(parts, SPLIT_BY_LENGTH_1)


@env.task
async def test_seqkit_split2_paired_end_length() -> None:
    parts = await seqkit_split2(await paired_reads(), args="--by-length 8K")
    expected = SPLIT_BY_LENGTH_1 | {
        "test_2.part_001.fastq.gz": "77b9076f82a762711582584342bde5a1",
        "test_2.part_002.fastq.gz": "33bb6e3edc759baa7ba6580da36def48",
    }
    await check_split2(parts, expected)


@env.task
async def test_seqkit_split2_single_end_part() -> None:
    parts = await seqkit_split2([await fixture(FASTQ + "test_1.fastq.gz")], args="--by-part 3")
    await check_split2(parts, SPLIT_BY_PART_1)


@env.task
async def test_seqkit_split2_paired_end_part() -> None:
    parts = await seqkit_split2(await paired_reads(), args="--by-part 3")
    expected = SPLIT_BY_PART_1 | {
        "test_2.part_001.fastq.gz": "13d14d69744bd93c8c51873b529bf714",
        "test_2.part_002.fastq.gz": "5f5a78d78f312164a1159c62d65c15f1",
        "test_2.part_003.fastq.gz": "48bcbceb485b73bcf1f198e252b016d8",
    }
    await check_split2(parts, expected)


@env.task
async def test_seqkit_split2_single_end_size() -> None:
    parts = await seqkit_split2([await fixture(FASTQ + "test_1.fastq.gz")], args="--by-size 50")
    await check_split2(parts, SPLIT_BY_SIZE_1)


@env.task
async def test_seqkit_split2_paired_end_size() -> None:
    parts = await seqkit_split2(await paired_reads(), args="--by-size 50")
    expected = SPLIT_BY_SIZE_1 | {
        "test_2.part_001.fastq.gz": "201ee95b559240e27830970b78a547c8",
        "test_2.part_002.fastq.gz": "35ff29a76f34b2507a37287352324650",
    }
    await check_split2(parts, expected)


@env.task
async def test_seqkit_split2_single_end_length_uncompressed() -> None:
    gzipped = await fixture(FASTQ + "test_1.fastq.gz")
    async with gzipped.open("rb") as fh:
        raw = bytes(await fh.read())
    path = Path(tempfile.mkdtemp()) / "test_1.fastq"
    path.write_bytes(gzip.decompress(raw))
    parts = await seqkit_split2([await File.from_local(str(path))], args="--by-length 8K")
    await check_split2(parts, {name.removesuffix(".gz"): md5 for name, md5 in SPLIT_BY_LENGTH_1.items()})


@env.task
async def test_seqkit_fq2fa_fastq_gz() -> None:
    out = await seqkit_fq2fa(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    assert name_of(out) == "test.fa.gz", f"expected test.fa.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, "f0c5c9110ce19e9ebbc9a6b6baf9e105", label="test.fa.gz")


HUMAN_GENOME = "genomics/homo_sapiens/genome/"
CHR22_MD5 = "f315020d899597c1b57e5fe9f60f4c3e"


async def chr22_header() -> File:
    path = Path(tempfile.mkdtemp()) / "chr22.header"
    path.write_text("chr22\n")
    return await File.from_local(str(path))


@env.task
async def test_seqkit_grep_file_gzipped() -> None:
    sequence = await fixture(HUMAN_GENOME + "genome.fasta.gz")
    out = await seqkit_grep(sequence, patterns=await chr22_header(), prefix="test")
    assert name_of(out) == "test.fa.gz", f"expected test.fa.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, CHR22_MD5, label="test.fa.gz")


@env.task
async def test_seqkit_grep_file_uncompressed() -> None:
    sequence = await fixture(HUMAN_GENOME + "genome.fasta")
    out = await seqkit_grep(sequence, patterns=await chr22_header(), prefix="test")
    assert name_of(out) == "test.fa", f"expected test.fa, got {name_of(out)}"
    await assert_md5(out, CHR22_MD5, label="test.fa")


@env.task
async def test_seqkit_grep_pattern() -> None:
    out = await seqkit_grep(await fixture(HUMAN_GENOME + "genome.fasta"), pattern="chr22", prefix="test")
    assert name_of(out) == "test.fa", f"expected test.fa, got {name_of(out)}"
    await assert_md5(out, CHR22_MD5, label="test.fa")


@env.task
async def test_seqkit_head_paired_fastq_gz() -> None:
    # upstream case: seqkit/head "sarscov2 - paired fastq.gz"
    reads = [await fixture(FASTQ + "test2_1.fastq.gz"), await fixture(FASTQ + "test2_2.fastq.gz")]
    subsets = await seqkit_head(reads, 100, prefix="test")
    expected = {
        "test_subset_test2_1.fastq.gz": "4161df271f9bfcd25d5845a1e220dbec",
        "test_subset_test2_2.fastq.gz": "2ebae722295ea66d84075a3b042e2b42",
    }
    assert [name_of(f) for f in subsets] == list(expected), f"wrote {[name_of(f) for f in subsets]}"
    for subset in subsets:
        await assert_gunzipped_md5(subset, expected[name_of(subset)], label=name_of(subset))


@env.task
async def test_seqkit_sliding_fasta() -> None:
    # upstream case: seqkit/sliding "sarscov2 - fasta"
    windows = await seqkit_sliding(await fixture(GENOME + "genome.fasta"), prefix="test", args="-s 2 -W 5")
    assert name_of(windows) == "test.fasta", f"wrote {name_of(windows)}"
    await assert_md5(windows, "d5e112b55c37bef1f16f28ba16b323dd", label="test.fasta")


@env.task
async def test_seqkit_sliding_fastq() -> None:
    # upstream case: seqkit/sliding "sarscov2 - fastq"
    windows = await seqkit_sliding(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", args="-s 2 -W 5")
    assert name_of(windows) == "test.fastq", f"wrote {name_of(windows)}"
    await assert_md5(windows, "59f475cac6d2c372ce2d3a1e3c40305b", label="test.fastq")


# upstream cases: seqkit/sort, with the test config's gzip ``-n``.
@env.task
async def test_seqkit_sort_genome_fasta() -> None:
    out = await seqkit_sort(await fixture(GENOME + "genome.fasta"), prefix="test", gzip_args="-n")
    assert name_of(out) == "test.fasta", f"expected test.fasta, got {name_of(out)}"
    await assert_md5(out, "483f4a5dfe60171c86ee9b7e6dff908b", label="test.fasta")


@env.task
async def test_seqkit_sort_genome_fasta_gz() -> None:
    out = await seqkit_sort(await fixture(GENOME + "genome.fasta.gz"), prefix="test", gzip_args="-n")
    assert name_of(out) == "test.fasta.gz", f"expected test.fasta.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, "483f4a5dfe60171c86ee9b7e6dff908b", label="test.fasta.gz")


@env.task
async def test_seqkit_sort_fastq_gz() -> None:
    out = await seqkit_sort(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test", gzip_args="-n")
    assert name_of(out) == "test.fastq.gz", f"expected test.fastq.gz, got {name_of(out)}"
    await assert_gunzipped_md5(out, "80d9fd1976648214f44e42bc3d705972", label="test.fastq.gz")


@env.task
async def test_seqkit_sort_name_conflict_fails() -> None:
    try:
        await seqkit_sort(await fixture(FASTQ + "test_1.fastq.gz"), prefix="test_1", gzip_args="-n")
    except ValueError as e:
        assert "Input and output names are the same" in str(e), str(e)
    else:
        raise AssertionError("seqkit_sort accepted an output name equal to its input's")


tests = [
    test_seqkit_stats_single_end,
    test_seqkit_stats_paired_end,
    test_seqkit_stats_nanopore,
    test_seqkit_stats_genome_fasta,
    test_seqkit_stats_transcriptome_fasta,
    test_seqkit_replace_fasta_gz,
    test_seqkit_replace_uncompressed,
    test_seqkit_replace_custom_extension,
    test_seqkit_seq_genome_fasta,
    test_seqkit_seq_genome_fasta_gz,
    test_seqkit_seq_fastq_gz,
    test_seqkit_seq_name_conflict_fails,
    test_seqkit_rmdup_genome_fasta,
    test_seqkit_rmdup_repeated_fasta,
    test_seqkit_rmdup_genome_fasta_gz,
    test_seqkit_rmdup_fastq_gz,
    test_seqkit_rmdup_name_conflict_fails,
    test_seqkit_split2_single_end_length,
    test_seqkit_split2_paired_end_length,
    test_seqkit_split2_single_end_part,
    test_seqkit_split2_paired_end_part,
    test_seqkit_split2_single_end_size,
    test_seqkit_split2_paired_end_size,
    test_seqkit_split2_single_end_length_uncompressed,
    test_seqkit_fq2fa_fastq_gz,
    test_seqkit_grep_file_gzipped,
    test_seqkit_grep_file_uncompressed,
    test_seqkit_grep_pattern,
    test_seqkit_head_paired_fastq_gz,
    test_seqkit_sliding_fasta,
    test_seqkit_sliding_fastq,
    test_seqkit_sort_genome_fasta,
    test_seqkit_sort_genome_fasta_gz,
    test_seqkit_sort_fastq_gz,
    test_seqkit_sort_name_conflict_fails,
]
