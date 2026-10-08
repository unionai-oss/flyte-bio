"""Tests for flyte_bio.modules.umitools.

Upstream cases from umitools/{extract,dedup,prepareforrsem}. Upstream's
extract snapshot records only the output names; prepareforrsem snapshots
the output BAM itself; dedup snapshots, with stats, the three TSVs and the
BAM's read sequences. The dedup BAM md5s come from the module's newer
snapshot (nf-core/modules 4dd270e, ``md5Reads``; same UMI-tools 1.1.6
container) rather than the one rnaseq pins, whose ``getSamLinesMD5``
hashes htsjdk's SAM text, which ``samtools view`` can't reproduce.
"""

from flyte_bio.modules.umitools import umitools_dedup, umitools_extract, umitools_prepareforrsem
from tests.framework import assert_md5, assert_nonempty, assert_reads_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
BAM = "genomics/sarscov2/illumina/bam/"


@env.task
async def test_extract_single_end() -> None:
    # upstream case: umitools/extract "single end" (ext.args --bc-pattern="NNNN")
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await umitools_extract(reads, prefix="test", args="--bc-pattern=NNNN")
    assert result.reads_2 is None
    assert result.reads_1.path.endswith("test.umi_extract.fastq.gz"), result.reads_1.path
    assert result.log.path.endswith("test.umi_extract.log"), result.log.path
    await assert_nonempty(result.reads_1, label="extracted reads")


@env.task
async def test_extract_paired_end() -> None:
    # upstream case: umitools/extract "pair end"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    result = await umitools_extract(r1, r2, prefix="test", args="--bc-pattern=NNNN")
    assert result.reads_2 is not None
    assert result.reads_1.path.endswith("test.umi_extract_1.fastq.gz"), result.reads_1.path
    assert result.reads_2.path.endswith("test.umi_extract_2.fastq.gz"), result.reads_2.path
    await assert_nonempty(result.reads_1, label="extracted R1")
    await assert_nonempty(result.reads_2, label="extracted R2")
    await assert_nonempty(result.log, label="extract log")


@env.task
async def test_dedup_single_end_no_stats() -> None:
    # upstream case: umitools/dedup "se - no stats" (ext.prefix "${meta.id}.dedup")
    bam = await fixture(BAM + "test.single_end.umi.sorted.bam")
    bai = await fixture(BAM + "test.single_end.umi.sorted.bam.bai")
    result = await umitools_dedup(bam, bai, prefix="test.dedup", paired=False)
    await assert_reads_md5(result.bam, "4ec386f9c27e2e9147f14f20fc0526e9", label="se dedup bam")
    await assert_nonempty(result.log, label="se dedup log")
    assert result.edit_distance is None


@env.task
async def test_dedup_paired_end_with_stats() -> None:
    # upstream cases: umitools/dedup "pe - no stats" and "pe - with stats" (same BAM)
    bam = await fixture(BAM + "test.paired_end.umi.sorted.bam")
    bai = await fixture(BAM + "test.paired_end.umi.sorted.bam.bai")
    result = await umitools_dedup(bam, bai, prefix="test.dedup", paired=True, output_stats=True)
    await assert_reads_md5(result.bam, "c1917631c47d16320d002b867e226a2e", label="pe dedup bam")
    assert result.edit_distance and result.per_umi and result.per_umi_per_position
    await assert_md5(result.edit_distance, "c247a49b58768e6e2e86a6c08483e612", label="edit_distance.tsv")
    await assert_md5(result.per_umi, "ced75f7bdbf38bf78f3137d5325a8773", label="per_umi.tsv")
    await assert_md5(result.per_umi_per_position, "2e1a12e6f720510880068deddeefe063", label="per_umi_per_position.tsv")


@env.task
async def test_prepareforrsem() -> None:
    # upstream case: umitools/prepareforrsem "sarscov2 - bam"
    bam = await fixture(BAM + "test.paired_end.sorted.bam")
    result = await umitools_prepareforrsem(bam, prefix="test")
    await assert_md5(result.bam, "fe4b8302182615651e4e7784ec67c819", label="prepare-for-rsem bam")
    await assert_nonempty(result.log, label="prepare-for-rsem log")


tests = [
    test_extract_single_end,
    test_extract_paired_end,
    test_dedup_single_end_no_stats,
    test_dedup_paired_end_with_stats,
    test_prepareforrsem,
]
