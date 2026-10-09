"""Tests for flyte_bio.modules.samtools.

md5 anchors come from the upstream samtools module snapshots (pinned to the
same tool version the RNA-seq pipeline uses). We md5 only pure-data outputs
— ``.fai``/``.sizes`` and the flagstat/idxstats count tables, which are
derived solely from input content. ``sort``/``index``/``stats`` produce
outputs that embed the command line, a version banner, or thread-dependent
compression, so they're checked run-to-green instead.
"""

from flyte_bio.modules.samtools import (
    samtools_faidx,
    samtools_flagstat,
    samtools_idxstats,
    samtools_index,
    samtools_merge,
    samtools_sort,
    samtools_stats,
    samtools_view,
)
from tests.framework import assert_md5, assert_nonempty, assert_reads_md5, env, fixture

SORTED_BAM = "genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"


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
async def test_index() -> None:
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


tests = [
    test_faidx,
    test_sort,
    test_view_primary_bam,
    test_merge_bams,
    test_index,
    test_stats,
    test_stats_cram,
    test_flagstat,
    test_idxstats,
]
