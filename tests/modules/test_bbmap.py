"""Tests for flyte_bio.modules.bbmap.

Mirrors the upstream bbmap/bbsplit "split with prebuilt index" case: build
an index, then split reads with it. Expected md5s are the upstream snapshot
values (nf-test hashes the decompressed reads). The log check proves the
staged index was reused rather than rebuilt — the timestamp/path fix-up
upstream applies to a staged index is what makes that work.

The bbmap/bbduk cases are upstream's (modules commit 9a1cc9f, BBMap 39.18),
with their test config's ``trimq=10 qtrim=r`` and ``<id>.trim`` prefix. The
contaminant cases' snapshots record only the version, so they check the
output names and the log, as upstream asserts.

The bbmap/bbnorm cases are upstream's at the same commit, with their test
config's ``target=100 min=5``. Upstream snapshots lines of the decompressed
output (all of it, the first two of each mate, or lines 3-7), so these tests
compare the same lines. Upstream gives BBNorm 95% of the memory as heap;
here it gets the usual 80%.

The bbmap/filterbyname cases are upstream's at the same commit. Each upstream
case declares ``single_end: false``, and although two of them list both
test_1.fastq.gz and test_2.fastq.gz, the process's ``tuple val(meta),
path(reads)`` input takes only the first file: every case filters
test_1.fastq.gz alone, read as interleaved pairs. These tests pass the same
single file, so the snapshot md5s reproduce (the ``.bz2`` outputs, written by
bzip2 -9, byte for byte). The paired-end bz2 case also runs a real pair of
mates, which upstream doesn't: with no names the reads pass unchanged.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.bbmap import (
    BBDukResult,
    BBNormResult,
    FilterByNameResult,
    bbduk,
    bbmap_filterbyname,
    bbnorm,
    bbsplit,
    bbsplit_index,
)
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
BBDUK_ARGS = "trimq=10 qtrim=r"


@env.task
async def test_bbsplit_prebuilt_index() -> None:
    # upstream case: bbmap/bbsplit "sarscov2_se_fastq_fasta_chr22_fasta - split with prebuilt index"
    primary = await fixture("genomics/sarscov2/genome/genome.fasta")
    human = await fixture("genomics/homo_sapiens/genome/chr22/sequence/chr22_23800000-23980000.fa")
    index = await bbsplit_index(primary, {"human": human})

    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    result = await bbsplit(reads, None, index, prefix="test")
    await assert_gunzipped_md5(result.reads_1, "4161df271f9bfcd25d5845a1e220dbec", label="primary reads")
    await assert_md5(result.stats, "2cbf69b72e5f4f8508306b54e8fe2861", label="stats")
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    assert "If you wish to regenerate the index" in log, "bbsplit rebuilt the index instead of reusing it"


async def check_bbduk(result: BBDukResult, names: list[str], md5s: list[str] | None, mode: str) -> None:
    assert [f.name for f in result.reads] == names, [f.name for f in result.reads]
    for file, md5 in zip(result.reads, md5s or [None] * len(names), strict=True):
        if md5 is None:
            await assert_nonempty(file, label=file.name)
        else:
            await assert_gunzipped_md5(file, md5, label=file.name)
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    assert f"Input is being processed as {mode}" in log, f"bbduk log lacks 'processed as {mode}'"


async def bbduk_reads() -> tuple[File, File]:
    return await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")


@env.task
async def test_bbduk_single_end() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - single end fastq - fastq"
    reads_1, _ = await bbduk_reads()
    result = await bbduk(reads_1, prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(result, ["test.trim.fastq.gz"], ["4161df271f9bfcd25d5845a1e220dbec"], "unpaired")


@env.task
async def test_bbduk_paired_end() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - paired end fastq - fastq"
    result = await bbduk(*await bbduk_reads(), prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(
        result,
        ["test.trim_1.fastq.gz", "test.trim_2.fastq.gz"],
        ["4161df271f9bfcd25d5845a1e220dbec", "2ebae722295ea66d84075a3b042e2b42"],
        "paired",
    )


@env.task
async def test_bbduk_single_end_contaminants() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - single end w/ contams [fastq,fasta] - fastq"
    reads_1, _ = await bbduk_reads()
    contaminants = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    result = await bbduk(reads_1, contaminants=contaminants, prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(result, ["test.trim.fastq.gz"], None, "unpaired")


@env.task
async def test_bbduk_paired_end_contaminants() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - paired end w/ contams [fastq,fasta] - fastq"
    contaminants = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    result = await bbduk(*await bbduk_reads(), contaminants=contaminants, prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(result, ["test.trim_1.fastq.gz", "test.trim_2.fastq.gz"], None, "paired")


BBNORM_ARGS = "target=100 min=5"


async def gunzipped_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode().splitlines()


async def check_bbnorm(result: BBNormResult, names: list[str]) -> None:
    assert [f.name for f in result.reads] == names, [f.name for f in result.reads]
    assert result.log.name == "test.bbnorm.log", result.log.name


@env.task
async def test_bbnorm_single_end() -> None:
    # upstream case: bbmap/bbnorm "test-bbmap-bbnorm-se" (every read falls below min depth)
    reads_1, _ = await bbduk_reads()
    result = await bbnorm(reads_1, prefix="test", args=BBNORM_ARGS)
    await check_bbnorm(result, ["test.fastq.gz"])
    lines = await gunzipped_lines(result.reads[0])
    assert lines == [], lines[:4]


@env.task
async def test_bbnorm_paired_end() -> None:
    # upstream case: bbmap/bbnorm "test-bbmap-bbnorm-pe"
    result = await bbnorm(*await bbduk_reads(), prefix="test", args=BBNORM_ARGS)
    await check_bbnorm(result, ["test_1.nm.fastq.gz", "test_2.nm.fastq.gz"])
    heads = [(await gunzipped_lines(f))[0:2] for f in result.reads]
    assert heads == [
        [
            "@ERR5069949.2161340 NS500628:121:HK3MMAFX2:2:21107:15810:18317/1",
            "AACCAGAATATTTCAATTCAGTGTGTAGACTTATGAAAACTATAGGTCCAGACATGTTCCTCGGAACTTGTCGGCGTTGT",
        ],
        [
            "@ERR5069949.2161340 NS500628:121:HK3MMAFX2:2:21107:15810:18317/2",
            "GGACAACGCCGACAAGTTCCGAGGAACATGTCTGGACCTATAGTTTTCATAAGTCTACACACTGAATTGAAATATTCTGGTT",
        ],
    ], heads


@env.task
async def test_bbnorm_interleaved() -> None:
    # upstream case: bbmap/bbnorm "test-bbmap-bbnorm-interleaved"
    reads = await fixture(FASTQ + "test_interleaved.fastq.gz")
    result = await bbnorm([reads], prefix="test", args=BBNORM_ARGS)
    await check_bbnorm(result, ["test.fastq.gz"])
    lines = (await gunzipped_lines(result.reads[0]))[3:8]
    assert lines == [
        "A/AA//EEAEA/E/AEEEE6EE/EEEA/6AEEEEEEEEE6EEEAEAEE//A/EEEEEE//E/E/A//E/E/<<EE</E/E",
        "@ERR5069949.2161340 NS500628:121:HK3MMAFX2:2:21107:15810:18317/2",
        "GGACAACGCCGACAAGTTCCGAGGAACATGTCTGGACCTATAGTTTTCATAAGTCTACACACTGAATTGAAATATTCTGGTT",
        "+",
        "AA6AAAE6/EA/<E/E/<//EE///EAEAEEEE////6EEE//EEE/EAEEE/AEA/EA<EA</AEEE/EE/AEAE</E//A",
    ], lines


@env.task
async def test_bbnorm_multiple_input() -> None:
    # upstream case: bbmap/bbnorm "test-bbmap-bbnorm-multiple-input" (two single-end FASTQs, read in turn)
    result = await bbnorm(list(await bbduk_reads()), prefix="test", args=BBNORM_ARGS)
    await check_bbnorm(result, ["test.fastq.gz"])
    lines = (await gunzipped_lines(result.reads[0]))[3:8]
    assert lines == [
        (
            "AAAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
            "EEEEEEEAAEEEEE<EEAAAEEEEEEEEEAAAAEAEEEAEEEEEE<AAAA"
        ),
        "@ERR5069949.1331889 NS500628:121:HK3MMAFX2:3:11605:17658:18222/1",
        (
            "GTCTACAAGCTGGTAATGCAACAGAAGTGCCTGCCAATTCAACTGTATTATCTTTCTGTGCTTTTGCTGTAGATGCTGCTAAAGCTTACAAAGATTATCTAGCTAGTG"
            "GGGGACAACCAATCACTAATTGTG"
        ),
        "+",
        (
            "A/AAAEEEEEEEEEEEEEEEEEEEAEEEEAEEEEEEEEEEEAEEEEEEEEEEEEEEE/EEEEE<AEAEEEEE/EAEAEEE/AEEEEEEEEEEEEEEEEEEEEAE/EEE"
            "EEEEEEEEEEEEEEEEEEEEA<EE"
        ),
    ], lines


FILTER_NAMES = "ERR5069949.2151832,ERR5069949.576388,ERR5069949.501486"


async def check_filterbyname(result: FilterByNameResult, md5s: dict[str, str], gunzipped: bool = False) -> None:
    assert [f.name for f in result.reads] == list(md5s), [f.name for f in result.reads]
    assert result.log.name == "test.log", result.log.name
    for file, md5 in zip(result.reads, md5s.values(), strict=True):
        if gunzipped:
            await assert_gunzipped_md5(file, md5, label=file.name)
        else:
            await assert_md5(file, md5, label=file.name)


@env.task
async def test_bbmap_filterbyname_paired_bz2() -> None:
    # upstream case: bbmap/filterbyname "paired end fastq.bz2" (test_1.fastq.gz read as interleaved pairs)
    reads_1, reads_2 = await bbduk_reads()
    result = await bbmap_filterbyname(reads_1, single_end=False, output_format="fastq.bz2", prefix="test")
    await check_filterbyname(
        result,
        {
            "test_1.fastq.bz2": "410105042b5c32ed304342f3259f7f80",
            "test_2.fastq.bz2": "8f7510b06e5123b9511f9ae54863f85d",
        },
    )
    # Not upstream: two mates, no names, gzipped output: the reads are unchanged.
    mates = await bbmap_filterbyname(reads_1, reads_2, prefix="test")
    await check_filterbyname(
        mates,
        {"test_1.fastq.gz": "4161df271f9bfcd25d5845a1e220dbec", "test_2.fastq.gz": "2ebae722295ea66d84075a3b042e2b42"},
        gunzipped=True,
    )


@env.task
async def test_bbmap_filterbyname_fasta() -> None:
    # upstream case: bbmap/filterbyname "single end fasta"
    reads_1, _ = await bbduk_reads()
    result = await bbmap_filterbyname(reads_1, single_end=False, output_format="fasta", prefix="test")
    await check_filterbyname(
        result, {"test_1.fasta": "884c105cedb985918d13266c8f2fc237", "test_2.fasta": "2c6b38825908c9bec9748f486d030487"}
    )


@env.task
async def test_bbmap_filterbyname_names() -> None:
    # upstream case: bbmap/filterbyname "single end fastq.gz filter"
    reads_1, _ = await bbduk_reads()
    result = await bbmap_filterbyname(
        reads_1, names=FILTER_NAMES, single_end=False, output_format="fasta", prefix="test"
    )
    await check_filterbyname(
        result, {"test_1.fasta": "8408578e9a56eb5c460304533544ea3f", "test_2.fasta": "582eff5fc721306af42bfd3718b61548"}
    )


@env.task
async def test_bbmap_filterbyname_interleaved_output() -> None:
    # upstream case: bbmap/filterbyname "paired end fastq.gz filter interleaved"
    reads_1, _ = await bbduk_reads()
    result = await bbmap_filterbyname(
        reads_1, names=FILTER_NAMES, single_end=False, interleaved_output=True, prefix="test"
    )
    await check_filterbyname(result, {"test.fastq.gz": "91dc9a4c9423b78c36bbce6b9c851b9b"}, gunzipped=True)


tests = [
    test_bbsplit_prebuilt_index,
    test_bbduk_single_end,
    test_bbduk_paired_end,
    test_bbduk_single_end_contaminants,
    test_bbduk_paired_end_contaminants,
    test_bbnorm_single_end,
    test_bbnorm_paired_end,
    test_bbnorm_interleaved,
    test_bbnorm_multiple_input,
    test_bbmap_filterbyname_paired_bz2,
    test_bbmap_filterbyname_fasta,
    test_bbmap_filterbyname_names,
    test_bbmap_filterbyname_interleaved_output,
]
