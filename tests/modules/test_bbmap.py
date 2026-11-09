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

The bbmap/repair cases are upstream's at the same commit, with their test
config's ``qin=33``. Upstream snapshots the repaired reads' md5s
(decompressed) and the singleton file's and log's names.

The bbmap/index, bbmap/align, bbmap/pileup and bbmap/clumpify cases are
upstream's at the same commit. bbmap/index snapshots the index's file names
and one stable file's md5. bbmap/align snapshots BAMs by their reads' md5,
and reads outputs by name, line count and the md5 of their sorted lines
(joined by newlines); its "pigz" case runs with ``unpigz=t``. bbmap/clumpify
snapshots the decompressed reads' md5s and the log's name.

The bbmap/sendsketch cases are upstream's at the same commit. They query JGI's
sketch server over the network, and upstream snapshots the md5 of its report.
The report reflects the server's live database: the paired case's no longer
matches upstream's snapshot (the hits' order and figures differ), and the
single case's will drift the same way. So both check the report's query line,
its columns and its hits instead.
"""

import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.bbmap import (
    BBDukResult,
    BBMapAlignResult,
    BBNormResult,
    FilterByNameResult,
    RepairResult,
    bbduk,
    bbmap_align,
    bbmap_clumpify,
    bbmap_filterbyname,
    bbmap_index,
    bbmap_pileup,
    bbmap_repair,
    bbmap_sendsketch,
    bbnorm,
    bbsplit,
    bbsplit_index,
)
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, assert_reads_md5, env, fixture

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


async def check_repair(result: RepairResult, md5s: dict[str, str]) -> None:
    assert [f.name for f in result.reads] == list(md5s), [f.name for f in result.reads]
    for file, md5 in zip(result.reads, md5s.values(), strict=True):
        await assert_gunzipped_md5(file, md5, label=file.name)
    assert result.singletons.name == "test_singleton.fastq.gz", result.singletons.name
    assert result.log.name == "test.repair.sh.log", result.log.name


@env.task
async def test_bbmap_repair_paired() -> None:
    # upstream case: bbmap/repair "sarscov2_illumina_paired - fastq_gz"
    result = await bbmap_repair(*await bbduk_reads(), prefix="test", args="qin=33")
    await check_repair(
        result,
        {
            "test_1_repaired.fastq.gz": "4161df271f9bfcd25d5845a1e220dbec",
            "test_2_repaired.fastq.gz": "2ebae722295ea66d84075a3b042e2b42",
        },
    )


@env.task
async def test_bbmap_repair_interleaved() -> None:
    # upstream case: bbmap/repair "sarscov2_illumina_interleaved - fastq_gz"
    reads = await fixture(FASTQ + "test_interleaved.fastq.gz")
    result = await bbmap_repair(reads, prefix="test", args="qin=33")
    await check_repair(result, {"test_repaired.fastq.gz": "efe40fe3df38b5a1c469890e5ad19874"})


HUMAN_GENOME = "genomics/homo_sapiens/genome/genome.fasta"
TRANSCRIPTOME = "genomics/sarscov2/genome/transcriptome.fasta"


@env.task
async def test_bbmap_index_human() -> None:
    # upstream case: bbmap/index "test-bbmap-index"
    index = await bbmap_index(await fixture(HUMAN_GENOME))
    files = [f async for f in index.walk()]
    names = sorted(f.name for f in files)
    assert names == [
        "chr1.chrom.gz",
        "chr1_index_k13_c15_b1.block",
        "chr1_index_k13_c15_b1.block2.gz",
        "info.txt",
        "scaffolds.txt.gz",
        "summary.txt",
    ], names
    block = next(f for f in files if f.name == "chr1_index_k13_c15_b1.block")
    await assert_md5(block, "9f0d9a7413c1d2c16cc24555b2381163", label=block.name)


async def check_align_bam(result: BBMapAlignResult, name: str, reads_md5: str) -> None:
    assert result.bam is not None and result.reads == [], "expected a BAM and no reads"
    assert result.bam.name == name, f"BAM named {result.bam.name}, expected {name}"
    assert result.log.name == "test.bbmap.log", result.log.name
    await assert_reads_md5(result.bam, reads_md5, label=name)


async def check_align_reads(result: BBMapAlignResult, expected: list[tuple[str, int, str]]) -> None:
    assert result.bam is None, "expected no BAM"
    assert [f.name for f in result.reads] == [name for name, _, _ in expected], [f.name for f in result.reads]
    assert result.log.name == "test.bbmap.log", result.log.name
    for file, (name, count, md5) in zip(result.reads, expected, strict=True):
        lines = await gunzipped_lines(file)
        actual = hashlib.md5("\n".join(sorted(lines)).encode()).hexdigest()
        assert (len(lines), actual) == (count, md5), f"{name}: {len(lines)} lines, sorted md5 {actual}"


PAIRED_BAM_READS_MD5 = "88d78879c34214288644aec7ec3cb270"


@env.task
async def test_bbmap_align_paired_fasta() -> None:
    # upstream case: bbmap/align "paired_end - fasta"
    result = await bbmap_align(*await bbduk_reads(), await fixture(HUMAN_GENOME), prefix="test")
    await check_align_bam(result, "test.bam", PAIRED_BAM_READS_MD5)


@env.task
async def test_bbmap_align_paired_index() -> None:
    # upstream case: bbmap/align "paired_end - index"
    index = await bbmap_index(await fixture(HUMAN_GENOME))
    result = await bbmap_align(*await bbduk_reads(), index, output_format="bam", output_reads="all", prefix="test")
    await check_align_bam(result, "test.bam", PAIRED_BAM_READS_MD5)


@env.task
async def test_bbmap_align_single_index() -> None:
    # upstream case: bbmap/align "single_end - index"
    index = await bbmap_index(await fixture(HUMAN_GENOME))
    reads_1, _ = await bbduk_reads()
    result = await bbmap_align(reads_1, None, index, output_format="bam", output_reads="all", prefix="test")
    await check_align_bam(result, "test.bam", "c17efa4ccc4e9018090792c71af92660")


@env.task
async def test_bbmap_align_paired_index_pigz() -> None:
    # upstream case: bbmap/align "paired_end - index - pigz"
    index = await bbmap_index(await fixture(HUMAN_GENOME))
    result = await bbmap_align(*await bbduk_reads(), index, prefix="test", args="unpigz=t")
    await check_align_bam(result, "test.bam", PAIRED_BAM_READS_MD5)


UNMAPPED_PAIR = [
    ("test.unmapped_1.fastq.gz", 4, "316e3bf71159aab93c75c735ab71c6ca"),
    ("test.unmapped_2.fastq.gz", 4, "e2c49ca4ab6f381c8f55c7a157e67e6a"),
]


@env.task
async def test_bbmap_align_transcriptome_all_fastq() -> None:
    # upstream case: bbmap/align "sarscov2 - [fastq1, fastq2] - transcriptome - all - fastq.gz"
    reads = await bbduk_reads()
    result = await bbmap_align(
        *reads, await fixture(TRANSCRIPTOME), output_format="fastq.gz", output_reads="all", prefix="test"
    )
    await check_align_reads(
        result,
        [
            ("test.all_1.fastq.gz", 400, "4de5088852e098ec6bd5dfb05425fe16"),
            ("test.all_2.fastq.gz", 400, "e79e5570171ab63573c2c704aa1c0c6b"),
        ],
    )


@env.task
async def test_bbmap_align_transcriptome_unmapped_fastq() -> None:
    # upstream case: bbmap/align "sarscov2 - [fastq1, fastq2] - transcriptome - unmapped - fastq.gz"
    reads = await bbduk_reads()
    result = await bbmap_align(
        *reads, await fixture(TRANSCRIPTOME), output_format="fastq.gz", output_reads="unmapped", prefix="test"
    )
    await check_align_reads(result, UNMAPPED_PAIR)


@env.task
async def test_bbmap_align_transcriptome_mapped_interleaved_output() -> None:
    # upstream case: bbmap/align "sarscov2 - [fastq1, fastq2] - transcriptome - mapped - fastq.gz - interleaved_output"
    reads = await bbduk_reads()
    result = await bbmap_align(
        *reads,
        await fixture(TRANSCRIPTOME),
        output_format="fastq.gz",
        output_reads="mapped",
        interleaved_output=True,
        prefix="test",
    )
    await check_align_reads(result, [("test.mapped.fastq.gz", 792, "48ab9757d7325cfa66b97e62b3d93073")])


@env.task
async def test_bbmap_align_interleaved_unmapped_fastq() -> None:
    # upstream case: bbmap/align "sarscov2 - interleaved fastq - transcriptome - unmapped - fastq.gz"
    reads = await fixture(FASTQ + "test_interleaved.fastq.gz")
    result = await bbmap_align(
        reads,
        None,
        await fixture(TRANSCRIPTOME),
        single_end=False,
        output_format="fastq.gz",
        output_reads="unmapped",
        prefix="test",
    )
    await check_align_reads(result, UNMAPPED_PAIR)


@env.task
async def test_bbmap_align_transcriptome_mapped_bam() -> None:
    # upstream case: bbmap/align "sarscov2 - fastq - transcriptome - mapped - bam"
    reads_1, _ = await bbduk_reads()
    result = await bbmap_align(
        reads_1, None, await fixture(TRANSCRIPTOME), output_format="bam", output_reads="mapped", prefix="test"
    )
    await check_align_bam(result, "test.mapped.bam", "946280c22bbf653e797187227a239974")


@env.task
async def test_bbmap_pileup_paired_bam() -> None:
    # upstream case: bbmap/pileup "test-bbmap-pileup"
    result = await bbmap_pileup(await fixture("genomics/sarscov2/illumina/bam/test.paired_end.bam"), prefix="test")
    assert result.covstats.name == "test.coverage.stats.txt", result.covstats.name
    assert result.hist.name == "test.coverage.hist.txt", result.hist.name
    await assert_md5(result.covstats, "c3fc9d0681589b69e3301ca3cb27b7a4", label=result.covstats.name)
    await assert_md5(result.hist, "96915920ef42ddc9483457dd4585a088", label=result.hist.name)


CLUMPED_1_MD5 = "27e51643262c1ef3905c4be184c3814c"


@env.task
async def test_bbmap_clumpify_single_end() -> None:
    # upstream case: bbmap/clumpify "test-bbmap-clumpify-single-end"
    reads_1, _ = await bbduk_reads()
    result = await bbmap_clumpify(reads_1, prefix="test")
    assert [f.name for f in result.reads] == ["test.clumped.fastq.gz"], [f.name for f in result.reads]
    assert result.log.name == "test.clumpify.log", result.log.name
    await assert_gunzipped_md5(result.reads[0], CLUMPED_1_MD5, label="test.clumped.fastq.gz")


@env.task
async def test_bbmap_clumpify_paired_end() -> None:
    # upstream case: bbmap/clumpify "test-bbmap-clumpify-paired-end"
    result = await bbmap_clumpify(*await bbduk_reads(), prefix="test")
    names = [f.name for f in result.reads]
    assert names == ["test_1.clumped.fastq.gz", "test_2.clumped.fastq.gz"], names
    assert result.log.name == "test.clumpify.log", result.log.name
    await assert_gunzipped_md5(result.reads[0], CLUMPED_1_MD5, label=names[0])
    await assert_gunzipped_md5(result.reads[1], "c70ab7bbd44d6b6fadd6a1a79ef1648f", label=names[1])


async def check_sendsketch_report(hits: File, query: str, expected_hit: str | None) -> None:
    """Check a sendsketch report's shape: the server's live database moves its hits and figures."""
    assert hits.name == "test.txt", f"wrote {hits.name}"
    async with hits.open("rb") as fh:
        report = bytes(await fh.read()).decode().splitlines()
    assert report[1].startswith(f"Query: {query}\tDB: RefSeq\t"), f"report's query line: {report[1]!r}"
    assert report[2].startswith("WKID\tKID\tANI\t"), f"report's header: {report[2]!r}"
    names = [line.split("\t")[-1] for line in report[3:] if line]
    assert names, "sendsketch reported no hits"
    if expected_hit is not None:
        assert expected_hit in names, f"no {expected_hit} hit among {names}"


@env.task
async def test_bbmap_sendsketch_single() -> None:
    # upstream case: bbmap/sendsketch "test-bbmap-sendsketch-single" (upstream's md5 ba6a7865... matched in 2026-11)
    hits = await bbmap_sendsketch(await fixture("genomics/homo_sapiens/genome/genome.fasta"), prefix="test")
    # The test genome is a small human slice; its top hits are other primates.
    await check_sendsketch_report(hits, "genome.fasta", expected_hit=None)


@env.task
async def test_bbmap_sendsketch_paired() -> None:
    # upstream case: bbmap/sendsketch "test-bbmap-sendsketch-paired": upstream sketches only the first mate
    hits = await bbmap_sendsketch(await fixture("genomics/homo_sapiens/illumina/fastq/test2_1.fastq.gz"), prefix="test")
    await check_sendsketch_report(hits, "test2_1.fastq.gz", expected_hit="Homo sapiens")


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
    test_bbmap_repair_paired,
    test_bbmap_repair_interleaved,
    test_bbmap_index_human,
    test_bbmap_align_paired_fasta,
    test_bbmap_align_paired_index,
    test_bbmap_align_single_index,
    test_bbmap_align_paired_index_pigz,
    test_bbmap_align_transcriptome_all_fastq,
    test_bbmap_align_transcriptome_unmapped_fastq,
    test_bbmap_align_transcriptome_mapped_interleaved_output,
    test_bbmap_align_interleaved_unmapped_fastq,
    test_bbmap_align_transcriptome_mapped_bam,
    test_bbmap_pileup_paired_bam,
    test_bbmap_clumpify_single_end,
    test_bbmap_clumpify_paired_end,
    test_bbmap_sendsketch_single,
    test_bbmap_sendsketch_paired,
]
