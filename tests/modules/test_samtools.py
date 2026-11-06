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
"""

import gzip

from flyte.io import File

from flyte_bio.modules.samtools import (
    samtools_collatefastq,
    samtools_faidx,
    samtools_fastq,
    samtools_flagstat,
    samtools_idxstats,
    samtools_index,
    samtools_merge,
    samtools_sort,
    samtools_stats,
    samtools_view,
)
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, assert_reads_md5, env, fixture

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
]
