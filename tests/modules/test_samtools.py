"""Tests for flyte_bio.modules.samtools.

md5 anchors come from the upstream samtools module snapshots (pinned to the
same tool version the RNA-seq pipeline uses). We md5 only pure-data outputs
— ``.fai``/``.sizes`` and the flagstat/idxstats count tables, which are
derived solely from input content. ``sort``/``index``/``stats`` produce
outputs that embed the command line, a version banner, or thread-dependent
compression, so they're checked run-to-green instead.

The FASTQ conversions' cases compare what upstream snapshots: the first
reads of each FASTQ (``samtools fastq``) and the decompressed FASTQs' md5s
(``samtools collatefastq``).

The convert and dict cases (upstream's samtools/convert and samtools/dict at
modules commit 9a1cc9f, samtools 1.25) check what upstream snapshots: the
output names, the converted BAM's read sequences' md5, and the dictionary's
first line.

The import cases are upstream's samtools/import at the same commit (samtools
1.25): the output name and the read sequences' md5, which is the same for a
pair and its interleaved file.

The trimheader case is upstream's samtools/trimheader at the same commit
(samtools 1.25): the trimmed header's @SQ lines, the read count and the read
sequences' md5.

The bam2fq and mpileup cases are upstream's at the same commit: the
decompressed FASTQs' and pileups' md5s (bam2fq's split case snapshots its
other and singleton files by name only).

The coverage, depth and reheader cases are upstream's at the same commit
(samtools 1.25), with their snapshots' md5s; coverage's two error cases
(a CRAM without its FASTA or its .fai) fail in the wrapper, as upstream's
module does before running.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.samtools import (
    samtools_bam2fq,
    samtools_calmd,
    samtools_collatefastq,
    samtools_convert,
    samtools_coverage,
    samtools_depth,
    samtools_dict,
    samtools_faidx,
    samtools_fastq,
    samtools_fixmate,
    samtools_flagstat,
    samtools_idxstats,
    samtools_import,
    samtools_index,
    samtools_markdup,
    samtools_merge,
    samtools_mpileup,
    samtools_reheader,
    samtools_sort,
    samtools_stats,
    samtools_trimheader,
    samtools_view,
)
from tests.framework import (
    assert_gunzipped_md5,
    assert_md5,
    assert_nonempty,
    assert_reads_md5,
    env,
    fixture,
    reads_md5,
)

SORTED_BAM = "genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"
PAIRED_BAM = "genomics/sarscov2/illumina/bam/test.paired_end.bam"
HUMAN_SORTED_BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"
HUMAN_SORTED_CRAM = "genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram"
HUMAN_FASTA = "genomics/homo_sapiens/genome/genome.fasta"

# The first 7 lines of `samtools fastq`'s R1 and R2 from PAIRED_BAM, as upstream snapshots them.
FASTQ_HEAD_1 = [
    "@ERR5069949.2151832",
    "TCATAAACCAAAGCACTCACAGTGTCAACAATTTCAGCAGGACAACGCCGACAAGTTCCGAGGAACATGTCTGGACCTATAGTTTTCATAAGTCTACACACTGAATTGAAATATTCTGGTTCTAGTGTGCCCTTAGTTAGCAATGTGCGT",
    "+",
    "AAAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEEEE<EEAAAEEEEEEEEEAAAAEAEEEAEEEEEE<AAAA",
    "@ERR5069949.576388",
    "ACTGTTTTCTTTGTAGAAAACATCCGTAATAGGACCTTTGTATTCTGAGGACTTTGTAAGTAAAGCACCGTCTATGC",
    "+",
]
FASTQ_HEAD_2 = [
    "@ERR5069949.2151832",
    "ATGTGTACATTGGCGACCCTGCTCAATTACCTGCACCACGCACATTGCTAACTAAGGGCACACTAGAACCAGAATATTTCAATTCAGTGTGTAGACTTATGAAAACTATAGGTCCAGACATGTTCCTCGGAACTTGTCGGCGTTGTCCTG",
    "+",
    "AAAAAEEEEEEEEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEEEEEEAEEEEEAAEEEEEEEEEAAEAAA<<EAAEEEEEEEAAA<<<AE",
    "@ERR5069949.576388",
    "GCATAGACGGTGCTTTACTTACAAAGTCCTCAGAATACAAAGGTCCTATTACGGATGTTTTCTACAAAGAAAACAGT",
    "+",
]
# The first 7 lines of the interleaved FASTQ.
INTERLEAVED_HEAD = [
    "@ERR5069949.2151832/1",
    FASTQ_HEAD_1[1],
    "+",
    FASTQ_HEAD_1[3],
    "@ERR5069949.2151832/2",
    FASTQ_HEAD_2[1],
    "+",
]


async def head_lines(file: File, count: int) -> list[str]:
    """The first ``count`` lines of ``file``, decompressed if it's gzipped."""
    async with file.open("rb") as fh:
        raw = bytes(await fh.read())
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode().splitlines()[:count]


def assert_name(file: File | None, name: str) -> None:
    assert file is not None and file.path.endswith("/" + name), f"expected {name}, got {file and file.path}"


@env.task
async def test_faidx() -> None:
    # upstream case: samtools/faidx "test_samtools_faidx_get_sizes"
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    fai, sizes = await samtools_faidx(fasta=fasta)
    await assert_md5(fai, "9da2a56e2853dc8c0b86a9e7229c9fe5", label="samtools faidx .fai")
    await assert_md5(sizes, "a57c401f27ae5133823fb09fb21c8a3c", label="samtools faidx .sizes")


@env.task
async def test_sort() -> None:
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.bam")
    out = await samtools_sort(bam=bam)
    await assert_nonempty(out, label="samtools sort")


@env.task
async def test_view_primary_bam() -> None:
    # The rnaseq UMI path's primary-only filter (upstream SAMTOOLS_VIEW_PRIMARY).
    bam = await fixture(SORTED_BAM)
    out = await samtools_view(bam=bam, args="-F 0x900 -b")
    async with out.open("rb") as fh:
        magic = bytes(await fh.read(2))
    assert magic == b"\x1f\x8b", f"samtools view -b wrote no BGZF BAM (starts {magic!r})"


@env.task
async def test_merge_bams() -> None:
    # upstream case: samtools/merge "bams" (md5 of the merged reads' sequences)
    bam = "genomics/sarscov2/illumina/bam/"
    inputs = [
        await fixture(bam + "test.paired_end.methylated.sorted.bam"),
        await fixture(bam + "test.paired_end.sorted.bam"),
        await fixture(bam + "test.single_end.sorted.bam"),
    ]
    merged = await samtools_merge(inputs, prefix="test")
    assert merged.path.endswith("/test.bam"), merged.path
    await assert_reads_md5(merged, "47c9f174d8c8afc1a13c75ee4b5e5d43", label="merged bams")


@env.task
async def test_samtools_index() -> None:
    bam = await fixture(SORTED_BAM)
    out = await samtools_index(bam=bam)
    await assert_nonempty(out, label="samtools index")


@env.task
async def test_stats() -> None:
    bam = await fixture(SORTED_BAM)
    out = await samtools_stats(bam=bam)
    await assert_nonempty(out, label="samtools stats")


@env.task
async def test_stats_cram() -> None:
    # upstream case: samtools/stats "cram" (its md5 covers a header with file paths, so check the content)
    hs = "genomics/homo_sapiens/"
    cram = await fixture(hs + "illumina/cram/test.paired_end.recalibrated.sorted.cram")
    fasta = await fixture(hs + "genome/chr21/sequence/genome.fasta")
    fai = await fixture(hs + "genome/chr21/sequence/genome.fasta.fai")
    out = await samtools_stats(bam=cram, fasta=fasta, fai=fai)
    async with out.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    assert text.startswith("# This file was produced by samtools stats"), text[:200]
    assert any(line.startswith("SN\traw total sequences:") for line in text.splitlines()), "no summary numbers"


@env.task
async def test_flagstat() -> None:
    # upstream case: samtools/flagstat "BAM" (note: the "BAM - stub" case's
    # md5 is a touched placeholder, not real output — don't use it)
    bam = await fixture(SORTED_BAM)
    out = await samtools_flagstat(bam=bam)
    await assert_md5(out, "4f7ffd1e6a5e85524d443209ac97d783", label="samtools flagstat")


@env.task
async def test_idxstats() -> None:
    # upstream case: samtools/idxstats "bam"
    bam = await fixture(SORTED_BAM)
    bai = await fixture(SORTED_BAM + ".bai")
    out = await samtools_idxstats(bam=bam, bai=bai)
    await assert_md5(out, "df60a8c8d6621100d05178c93fb053a2", label="samtools idxstats")


@env.task
async def test_samtools_fastq_bam() -> None:
    # upstream case: samtools/fastq "bam"
    out = await samtools_fastq(await fixture(PAIRED_BAM), prefix="test")
    assert len(out.fastq) == 2, out.fastq
    assert_name(out.fastq[0], "test_1.fastq.gz")
    assert_name(out.fastq[1], "test_2.fastq.gz")
    assert await head_lines(out.fastq[0], 7) == FASTQ_HEAD_1, "R1 reads differ from upstream's"
    assert await head_lines(out.fastq[1], 7) == FASTQ_HEAD_2, "R2 reads differ from upstream's"
    assert out.interleaved is None
    assert_name(out.singleton, "test_singleton.fastq.gz")
    assert_name(out.other, "test_other.fastq.gz")


@env.task
async def test_samtools_fastq_bam_interleave() -> None:
    # upstream case: samtools/fastq "bam_interleave"
    out = await samtools_fastq(await fixture(PAIRED_BAM), prefix="test", interleave=True)
    assert out.fastq == [] and out.singleton is None
    assert_name(out.interleaved, "test_interleaved.fastq")
    assert await head_lines(out.interleaved, 7) == INTERLEAVED_HEAD, "interleaved reads differ from upstream's"
    assert_name(out.other, "test_other.fastq.gz")


@env.task
async def test_samtools_collatefastq_bam_paired_end() -> None:
    # upstream case: samtools/collatefastq "human - bam - paired_end"
    out = await samtools_collatefastq(await fixture(HUMAN_SORTED_BAM), prefix="test", fasta=await fixture(HUMAN_FASTA))
    assert_name(out.fastq[0], "test_1.fq.gz")
    assert_name(out.fastq[1], "test_2.fq.gz")
    await assert_gunzipped_md5(out.fastq[0], "1cf671980643af6c1148ae5e8e94e350", label="collatefastq R1")
    await assert_gunzipped_md5(out.fastq[1], "38c1e9829115f9025f95435c5a4373d3", label="collatefastq R2")


@env.task
async def test_samtools_collatefastq_bam_single_end() -> None:
    # upstream case: samtools/collatefastq "human - bam - single_end" (snapshots only the "other" reads)
    out = await samtools_collatefastq(
        await fixture("genomics/sarscov2/illumina/bam/test.single_end.bam"),
        prefix="test",
        fasta=await fixture("genomics/sarscov2/genome/genome.fasta"),
        single_end=True,
    )
    assert len(out.fastq) == 1, out.fastq
    assert_name(out.other, "test_other.fq.gz")
    await assert_gunzipped_md5(out.other, "a6c101a06b5c9d5f8b91c0acd4ac5045", label="collatefastq other")


@env.task
async def test_samtools_collatefastq_cram() -> None:
    # upstream case: samtools/collatefastq "human - cram"
    out = await samtools_collatefastq(await fixture(HUMAN_SORTED_CRAM), prefix="test", fasta=await fixture(HUMAN_FASTA))
    await assert_gunzipped_md5(out.fastq[0], "1cf671980643af6c1148ae5e8e94e350", label="collatefastq cram R1")
    await assert_gunzipped_md5(out.fastq[1], "38c1e9829115f9025f95435c5a4373d3", label="collatefastq cram R2")


@env.task
async def test_samtools_collatefastq_bam_interleaved() -> None:
    # upstream case: samtools/collatefastq "human - bam - paired_end - interleaved"
    out = await samtools_collatefastq(
        await fixture(HUMAN_SORTED_BAM), prefix="test", fasta=await fixture(HUMAN_FASTA), interleave=True
    )
    assert out.fastq == [] and out.singleton is None
    assert_name(out.interleaved, "test_interleaved.fq")
    await assert_md5(out.interleaved, "4f2b93d492f0442fa89b02532c9b3530", label="collatefastq interleaved")


async def cram_reads_md5(cram: File, fasta: File, fai: File) -> str:
    """:func:`tests.framework.reads_md5` of a CRAM, decoded against its reference."""
    sam = await samtools_view(bam=cram, fasta=fasta, fai=fai)
    h = hashlib.md5()
    async with sam.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    for line in text.splitlines():
        h.update(line.split("\t")[9].encode())
    return h.hexdigest().lstrip("0")


@env.task
async def test_samtools_convert_bam_to_cram() -> None:
    # upstream case: samtools/convert "sarscov2 - [bam, bai], [fasta, fai]" (names only); that the CRAM
    # holds the BAM's reads is our check
    bam, bai, fasta, fai = await asyncio.gather(
        fixture(SORTED_BAM),
        fixture(SORTED_BAM + ".bai"),
        fixture("genomics/sarscov2/genome/genome.fasta"),
        fixture("genomics/sarscov2/genome/genome.fasta.fai"),
    )
    out = await samtools_convert(bam, fasta, "test", index=bai, fai=fai)
    assert_name(out.alignment, "test.cram")
    assert_name(out.index, "test.cram.crai")
    expected = await reads_md5(bam)
    actual = await cram_reads_md5(out.alignment, fasta, fai)
    assert actual == expected, f"CRAM reads md5 {actual}, expected the BAM's {expected}"


@env.task
async def test_samtools_convert_cram_to_bam() -> None:
    # upstream case: samtools/convert "homo_sapiens - [cram, crai], [fasta, fai]"
    cram, crai, fasta, fai = await asyncio.gather(
        fixture(HUMAN_SORTED_CRAM),
        fixture(HUMAN_SORTED_CRAM + ".crai"),
        fixture(HUMAN_FASTA),
        fixture(HUMAN_FASTA + ".fai"),
    )
    out = await samtools_convert(cram, fasta, "test", index=crai, fai=fai)
    assert_name(out.alignment, "test.bam")
    assert_name(out.index, "test.bam.bai")
    await assert_reads_md5(out.alignment, "2f11e4fe3390b8ad0a1852616fd1da04", label="converted bam")


@env.task
async def test_samtools_dict_fasta() -> None:
    # upstream case: samtools/dict "sarscov2 - fasta" (the first line); the @SQ line's name, length and
    # md5 are our check (its UR records the staged FASTA's path)
    out = await samtools_dict(await fixture("genomics/sarscov2/genome/genome.fasta"))
    assert_name(out, "genome.fasta.dict")
    lines = await head_lines(out, 2)
    assert lines[0] == "@HD\tVN:1.0\tSO:unsorted", lines[0]
    assert lines[1].startswith("@SQ\tSN:MT192765.1\tLN:29829\tM5:c95f3e5592d0ad9974e41e7f0ea14eb0\t"), lines[1]


async def sam_text(alignment: File, args: str) -> str:
    """``samtools view`` of ``alignment`` with ``args``, as text."""
    out = await samtools_view(bam=alignment, args=args)
    async with out.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_samtools_calmd_bam() -> None:
    # upstream case: samtools/calmd "bam" (modules commit 9a1cc9f, samtools 1.25). Its md5 is of calmd's SAM
    # text, @PG command line included: the 4 CPUs (-@ 4) and the inputs' names match upstream's.
    bam, fasta, fai = await asyncio.gather(
        fixture(HUMAN_SORTED_BAM), fixture(HUMAN_FASTA), fixture(HUMAN_FASTA + ".fai")
    )
    out = await samtools_calmd(bam, fasta, "test", fai=fai)
    assert_name(out, "test.bam")
    await assert_md5(out, "45f5722212b78c74427ac14670c7ea0c", label="samtools calmd")


@env.task
async def test_samtools_fixmate_bam() -> None:
    # upstream case: samtools/fixmate "sarscov2 - bam" (modules commit 9a1cc9f; upstream passes the human FASTA)
    bam, fasta, fai = await asyncio.gather(fixture(PAIRED_BAM), fixture(HUMAN_FASTA), fixture(HUMAN_FASTA + ".fai"))
    out = await samtools_fixmate(bam, "test_fixmate", fasta=fasta, fai=fai)
    assert_name(out, "test_fixmate.bam")
    await assert_reads_md5(out, "762e859a3d0ed1553655cde77665c940", label="samtools fixmate")


@env.task
async def test_samtools_fixmate_cram() -> None:
    # upstream case: samtools/fixmate "human - cram" (the name only). Upstream name-sorts the CRAM to a CRAM
    # first; our sort takes no reference, so the CRAM is converted to BAM and name-sorted. That the CRAM
    # holds the sorted reads is our check.
    cram, crai, fasta, fai = await asyncio.gather(
        fixture(HUMAN_SORTED_CRAM),
        fixture(HUMAN_SORTED_CRAM + ".crai"),
        fixture(HUMAN_FASTA),
        fixture(HUMAN_FASTA + ".fai"),
    )
    converted = await samtools_convert(cram, fasta, "test", index=crai, fai=fai)
    name_sorted = await samtools_sort(bam=converted.alignment, args="-n")
    out = await samtools_fixmate(
        name_sorted, "test_fixmate", fasta=fasta, fai=fai, args="--output-fmt cram,version=3.0"
    )
    assert_name(out, "test_fixmate.cram")
    expected = await reads_md5(name_sorted)
    actual = await cram_reads_md5(out, fasta, fai)
    assert actual == expected, f"fixmate CRAM reads md5 {actual}, expected the sorted input's {expected}"


@env.task
async def test_samtools_markdup_bam() -> None:
    # upstream case: samtools/markdup "sarscov2 - bam" (modules commit 9a1cc9f). Its md5 is of the whole BAM,
    # whose header records each step's command line, so check the content: the same reads in the same
    # order, and markdup's @PG line. (This fixture's reads have no duplicates to mark.)
    bam, fasta, fai = await asyncio.gather(
        fixture(PAIRED_BAM),
        fixture("genomics/sarscov2/genome/genome.fasta"),
        fixture("genomics/sarscov2/genome/genome.fasta.fai"),
    )
    fixed = await samtools_fixmate(bam, "test.fixmate", fasta=fasta, fai=fai)
    coordinate_sorted = await samtools_sort(bam=fixed)
    out = await samtools_markdup(coordinate_sorted, "test.markdup", fasta=fasta, fai=fai)
    assert_name(out, "test.markdup.bam")
    expected = await reads_md5(coordinate_sorted)
    await assert_reads_md5(out, expected, label="samtools markdup")
    header = await sam_text(out, "-H")
    assert "\tCL:samtools markdup " in header, header


SC2_FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_samtools_import_single() -> None:
    # upstream case: samtools/import "samtools_import_single "
    out = await samtools_import([await fixture(SC2_FASTQ + "test_1.fastq.gz")], "test", single_end=True)
    assert_name(out, "test.bam")
    await assert_reads_md5(out, "c17efa4ccc4e9018090792c71af92660", label="samtools import single-end")


@env.task
async def test_samtools_import_paired() -> None:
    # upstream case: samtools/import "samtools_import_paired "
    reads = list(await asyncio.gather(fixture(SC2_FASTQ + "test_1.fastq.gz"), fixture(SC2_FASTQ + "test_2.fastq.gz")))
    out = await samtools_import(reads, "test")
    assert_name(out, "test.bam")
    await assert_reads_md5(out, "88d78879c34214288644aec7ec3cb270", label="samtools import paired-end")


@env.task
async def test_samtools_import_interleaved() -> None:
    # upstream case: samtools/import "samtools_import_interleaved"
    out = await samtools_import([await fixture(SC2_FASTQ + "test_interleaved.fastq.gz")], "test")
    assert_name(out, "test.bam")
    await assert_reads_md5(out, "88d78879c34214288644aec7ec3cb270", label="samtools import interleaved")


@env.task
async def test_samtools_trimheader_bam() -> None:
    # upstream case: samtools/trimheader "bacteroides_fragilis - bam". The plasmid
    # NZ_CP069564.1 has no reads, and 47038 reads are unplaced.
    bam_path = "genomics/prokaryotes/bacteroides_fragilis/illumina/bam/test1.sorted.bam"
    out = await samtools_trimheader(await fixture(bam_path), await fixture(bam_path + ".bai"), "test")
    assert_name(out, "test.bam")
    header = await samtools_view(bam=out, args="-H")
    async with header.open("rb") as fh:
        sq = [line for line in bytes(await fh.read()).decode().splitlines() if line.startswith("@SQ")]
    assert sq == ["@SQ\tSN:NZ_CP069563.1\tLN:5150703"], f"trimmed @SQ lines: {sq}"
    count = await samtools_view(bam=out, args="-c")
    async with count.open("rb") as fh:
        reads = int(bytes(await fh.read()).decode().strip())
    assert reads == 100000, f"{reads} reads, expected 100000"
    await assert_reads_md5(out, "f18cb59ef156ae954b5939916b36f524", label="samtools trimheader")


UMI_CONVERTED_BAM = "genomics/homo_sapiens/illumina/bam/umi/test.paired_end.umi_converted.bam"


@env.task
async def test_samtools_bam2fq_bam() -> None:
    # upstream case: samtools/bam2fq "bam" (ext.args "-T RX")
    out = await samtools_bam2fq(await fixture(UMI_CONVERTED_BAM), prefix="test", args="-T RX")
    assert [f.name for f in out.fastq] == ["test_interleaved.fq.gz"], out.fastq
    assert out.other is None and out.singleton is None
    await assert_gunzipped_md5(out.fastq[0], "d733e66d29a4b366bf9df8c42f845256", label="bam2fq interleaved")


@env.task
async def test_samtools_bam2fq_bam_split() -> None:
    # upstream case: samtools/bam2fq "bam_split" (ext.args "-T RX")
    out = await samtools_bam2fq(await fixture(UMI_CONVERTED_BAM), prefix="test", split=True, args="-T RX")
    assert [f.name for f in out.fastq] == ["test_1.fq.gz", "test_2.fq.gz"], out.fastq
    assert_name(out.other, "test_other.fq.gz")
    assert_name(out.singleton, "test_singleton.fq.gz")
    await assert_gunzipped_md5(out.fastq[0], "2e4f2a6f06308e4b3eb8004b4044c949", label="bam2fq R1")
    await assert_gunzipped_md5(out.fastq[1], "bd81294b4907edab904d9692d79096f1", label="bam2fq R2")


@env.task
async def test_samtools_mpileup() -> None:
    # upstream cases: samtools/mpileup "bam - no fasta", "bam - fasta" and "bam - region - bai"
    sc2 = "genomics/sarscov2/"
    bam, bai, fasta, fai = await asyncio.gather(
        fixture(SORTED_BAM),
        fixture(SORTED_BAM + ".bai"),
        fixture(sc2 + "genome/genome.fasta"),
        fixture(sc2 + "genome/genome.fasta.fai"),
    )
    no_fasta, with_fasta, region = await asyncio.gather(
        samtools_mpileup(bam, prefix="test"),
        samtools_mpileup(bam, prefix="test", fasta=fasta, fai=fai),
        samtools_mpileup(bam, prefix="test", index=bai, args="-r MT192765.1"),
    )
    for out, label, md5 in (
        (no_fasta, "no fasta", "72340047c1feda82fa632d1762080bd8"),
        (with_fasta, "fasta", "958e6bead4103d72026f80153b6b5150"),
        (region, "region - bai", "72340047c1feda82fa632d1762080bd8"),
    ):
        assert_name(out, "test.mpileup.gz")
        await assert_gunzipped_md5(out, md5, label=f"mpileup {label}")


@env.task
async def test_samtools_coverage() -> None:
    # upstream cases: samtools/coverage "test_samtools_coverage_bam" and "test_samtools_coverage_cram"
    bam, bai, cram, crai, fasta, fai = await asyncio.gather(
        fixture(SORTED_BAM),
        fixture(SORTED_BAM + ".bai"),
        fixture(HUMAN_SORTED_CRAM),
        fixture(HUMAN_SORTED_CRAM + ".crai"),
        fixture(HUMAN_FASTA),
        fixture(HUMAN_FASTA + ".fai"),
    )
    from_bam, from_cram = await asyncio.gather(
        samtools_coverage(bam, "test", index=bai),
        samtools_coverage(cram, "test", index=crai, fasta=fasta, fai=fai),
    )
    assert_name(from_bam, "test.txt")
    await assert_md5(from_bam, "99a521b3bf53b6acf8055a44a571ea84", label="coverage bam")
    await assert_md5(from_cram, "ce896534bac51cfcc97e5508ae907e99", label="coverage cram")


@env.task
async def test_samtools_coverage_cram_without_reference() -> None:
    # upstream cases: samtools/coverage "test_samtools_coverage_cram_no_fasta" and "..._no_fai"
    cram, crai, fasta, fai = await asyncio.gather(
        fixture(HUMAN_SORTED_CRAM),
        fixture(HUMAN_SORTED_CRAM + ".crai"),
        fixture("genomics/sarscov2/genome/genome.fasta"),
        fixture("genomics/sarscov2/genome/genome.fasta.fai"),
    )
    for label, reference in (("no fasta", {"fai": fai}), ("no fai", {"fasta": fasta})):
        try:
            await samtools_coverage(cram, "test", index=crai, **reference)
        except ValueError as e:
            assert "no reference FASTA" in str(e), f"{label}: unexpected error {e}"
        else:
            raise AssertionError(f"samtools coverage accepted a CRAM with {label}")


@env.task
async def test_samtools_depth() -> None:
    # upstream cases: samtools/depth "test-samtools-depth" and "test-samtools-depth-intervals"
    bam, intervals = await asyncio.gather(
        fixture("genomics/homo_sapiens/illumina/bam/test.paired_end.recalibrated.sorted.bam"),
        fixture("genomics/homo_sapiens/genome/chr21/sequence/multi_intervals.bed"),
    )
    everywhere, at_intervals = await asyncio.gather(
        samtools_depth(bam, "test"), samtools_depth(bam, "test", intervals=intervals)
    )
    assert_name(everywhere, "test.tsv")
    await assert_md5(everywhere, "722d2c060944d6b96ae5ae00d8074657", label="depth")
    await assert_md5(at_intervals, "a89e2f8e1b66356c7500599fcadb8784", label="depth intervals")


@env.task
async def test_samtools_reheader_bam() -> None:
    # upstream cases: samtools/reheader "homo_sapiens - chrdel" and "homo_sapiens - rgdel"
    bam = await fixture(HUMAN_SORTED_BAM)
    chrdel, rgdel = await asyncio.gather(
        samtools_reheader(bam, "test_reheaded", command=r'sed -E "s/^(@SQ.*)(\tSN:)chr/\1\2/"'),
        samtools_reheader(bam, "test_reheaded", command="grep -v ^@RG"),
    )
    assert_name(chrdel, "test_reheaded.bam")
    await assert_md5(chrdel, "a89a52d5d858d5f0c9c7881dcf9cf20d", label="reheader chrdel")
    await assert_md5(rgdel, "e0276cbfc32ab6acebe7f3a3c03da7ae", label="reheader rgdel")


tests = [
    test_faidx,
    test_sort,
    test_view_primary_bam,
    test_merge_bams,
    test_samtools_index,
    test_stats,
    test_stats_cram,
    test_flagstat,
    test_idxstats,
    test_samtools_fastq_bam,
    test_samtools_fastq_bam_interleave,
    test_samtools_collatefastq_bam_paired_end,
    test_samtools_collatefastq_bam_single_end,
    test_samtools_collatefastq_cram,
    test_samtools_collatefastq_bam_interleaved,
    test_samtools_convert_bam_to_cram,
    test_samtools_convert_cram_to_bam,
    test_samtools_dict_fasta,
    test_samtools_calmd_bam,
    test_samtools_fixmate_bam,
    test_samtools_fixmate_cram,
    test_samtools_markdup_bam,
    test_samtools_import_single,
    test_samtools_import_paired,
    test_samtools_import_interleaved,
    test_samtools_trimheader_bam,
    test_samtools_bam2fq_bam,
    test_samtools_bam2fq_bam_split,
    test_samtools_mpileup,
    test_samtools_coverage,
    test_samtools_coverage_cram_without_reference,
    test_samtools_depth,
    test_samtools_reheader_bam,
]
