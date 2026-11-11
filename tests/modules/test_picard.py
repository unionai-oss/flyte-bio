"""Tests for flyte_bio.modules.picard.

Upstream's picard/markduplicates snapshot only records the output name and
the metrics header (which embeds run-specific paths), so this checks the
same shape: the marked BAM is produced and the metrics file is a Picard
duplication-metrics report for this command.

Upstream's picard/collectmultiplemetrics snapshots (modules commit 9a1cc9f,
Picard 3.5.0) record only the sorted names of the metrics files and PDFs,
since every metrics file's header records the command line and date. These
tests check the same names, plus that each metrics file is a Picard report of
this command.

The picard/addorreplacereadgroups, picard/createsequencedictionary and
picard/collecthsmetrics cases are upstream's at the same commit (Picard
3.5.0):

- AddOrReplaceReadGroups: the reads md5 of the output BAM/CRAM, with the
  test configs' args and ``<id>.replaced`` prefix, and the ``.bai`` name.
- CreateSequenceDictionary: upstream snapshots only the name, since the
  dictionary records the FASTA's absolute path (``UR:``). These tests also
  check its records against the test data's ``genome.dict``, ignoring ``UR:``.
- CollectHsMetrics: upstream snapshots the file's size and its first 101
  lines less the date line; these tests compare the md5 of those 100 lines
  (joined by newlines), and the size.
- MergeSamFiles (picard/mergesamfiles at the same commit): the merged BAM's
  name and reads md5.
- SortSam (picard/sortsam at the same commit): the sorted BAM's name and
  reads md5, with the test config's ``<id>.sorted`` prefix.
- CollectWgsMetrics (picard/collectwgsmetrics at the same commit, its own
  Picard 3.5.0 + R image): upstream snapshots only whether the metrics contain
  ``coverage    high_quality_coverage_count`` (spaces), which is ``false``, as
  the table is tab-separated. These tests check that, plus that the file is a
  Picard WGS-metrics report of this command with the tab-separated histogram.
- CollectInsertSizeMetrics: the metrics file's lines 6-9 (the metrics table)
  and the histogram's name.
- CollectRnaSeqMetrics: that the metrics name CollectRnaSeqMetrics and that no
  chart is written, with the refFlat from ucsc/gtftogenepred and the test
  config's args, as upstream.
- BedToIntervalList: the interval list's md5.
- LiftoverVcf (picard/liftovervcf at the same commit, with the test config's
  ``--WARN_ON_MISSING_CONTIG true``): nft-vcf's summary of the lifted VCF (no
  records) and ``variantsMD5`` of the unlifted one.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.picard import (
    CollectMultipleMetricsResult,
    picard_addorreplacereadgroups,
    picard_bedtointervallist,
    picard_collecthsmetrics,
    picard_collectinsertsizemetrics,
    picard_collectmultiplemetrics,
    picard_collectrnaseqmetrics,
    picard_collectwgsmetrics,
    picard_createsequencedictionary,
    picard_filtersamreads,
    picard_liftovervcf,
    picard_markduplicates,
    picard_mergesamfiles,
    picard_sortsam,
)
from flyte_bio.modules.samtools import samtools_view
from flyte_bio.modules.ucsc import gtftogenepred
from tests.framework import assert_md5, assert_nonempty, assert_reads_md5, env, fixture

SC2 = "genomics/sarscov2/"
COLLECTMULTIPLEMETRICS_METRICS = [
    "test.CollectMultipleMetrics.alignment_summary_metrics",
    "test.CollectMultipleMetrics.base_distribution_by_cycle_metrics",
    "test.CollectMultipleMetrics.insert_size_metrics",
    "test.CollectMultipleMetrics.quality_by_cycle_metrics",
    "test.CollectMultipleMetrics.quality_distribution_metrics",
]
COLLECTMULTIPLEMETRICS_PDFS = [
    "test.CollectMultipleMetrics.base_distribution_by_cycle.pdf",
    "test.CollectMultipleMetrics.insert_size_histogram.pdf",
    "test.CollectMultipleMetrics.quality_by_cycle.pdf",
    "test.CollectMultipleMetrics.quality_distribution.pdf",
    "test.CollectMultipleMetrics.read_length_histogram.pdf",
]

COMMAND = "# CollectMultipleMetrics --INPUT test.paired_end.sorted.bam "


async def text_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


@env.task
async def test_markduplicates_sorted_bam() -> None:
    # upstream case: picard/markduplicates "sarscov2 [sorted bam]"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam")
    result = await picard_markduplicates(bam, prefix="test.md", args="--ASSUME_SORT_ORDER queryname")
    await assert_nonempty(result.bam, label="marked bam")
    lines = await text_lines(result.metrics)
    assert lines[0] == "## htsjdk.samtools.metrics.StringHeader", lines[:2]
    assert lines[1].startswith("# MarkDuplicates --INPUT ") and "--OUTPUT " in lines[1], lines[1]
    assert any(line.startswith("## METRICS CLASS\tpicard.sam.DuplicationMetrics") for line in lines), "no metrics table"


async def check_collectmultiplemetrics(result: CollectMultipleMetricsResult, reference: bool) -> None:
    assert [f.name for f in result.metrics] == COLLECTMULTIPLEMETRICS_METRICS, [f.name for f in result.metrics]
    assert [f.name for f in result.pdfs] == COLLECTMULTIPLEMETRICS_PDFS, [f.name for f in result.pdfs]
    for metrics in result.metrics:
        lines = await text_lines(metrics)
        assert lines[0] == "## htsjdk.samtools.metrics.StringHeader", (metrics.name, lines[:2])
        command = lines[1]
        assert command.startswith(COMMAND), (metrics.name, command)
        assert ("--REFERENCE_SEQUENCE genome.fasta" in command) == reference, (metrics.name, command)
    for pdf in result.pdfs:
        async with pdf.open("rb") as fh:
            assert bytes(await fh.read(5)) == b"%PDF-", f"{pdf.name} is not a PDF"


async def sorted_bam() -> tuple[File, File]:
    bam = await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam")
    return bam, await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai")


@env.task
async def test_picard_collectmultiplemetrics() -> None:
    # upstream case: picard/collectmultiplemetrics "test-picard-collectmultiplemetrics"
    bam, bai = await sorted_bam()
    fasta = await fixture(SC2 + "genome/genome.fasta")
    result = await picard_collectmultiplemetrics(bam, prefix="test", bai=bai, fasta=fasta)
    await check_collectmultiplemetrics(result, reference=True)


@env.task
async def test_picard_collectmultiplemetrics_nofasta() -> None:
    # upstream case: picard/collectmultiplemetrics "test-picard-collectmultiplemetrics-nofasta"
    bam, bai = await sorted_bam()
    result = await picard_collectmultiplemetrics(bam, prefix="test", bai=bai)
    await check_collectmultiplemetrics(result, reference=False)


@env.task
async def test_picard_collectmultiplemetrics_cram() -> None:
    # upstream case: picard/collectmultiplemetrics "test-picard-collectmultiplemetrics-cram". Despite its
    # name, upstream passes the BAM, with the sarscov2 FASTA and a homo_sapiens .fai beside it.
    bam, bai = await sorted_bam()
    fasta = await fixture(SC2 + "genome/genome.fasta")
    fai = await fixture("genomics/homo_sapiens/genome/genome.fasta.fai")
    result = await picard_collectmultiplemetrics(bam, prefix="test", bai=bai, fasta=fasta, fai=fai)
    await check_collectmultiplemetrics(result, reference=True)


READ_GROUP_ARGS = "-LB test -PL ILLUMINA -PU bc1 -SM test"
HUMAN = "genomics/homo_sapiens/"


@env.task
async def test_picard_addorreplacereadgroups_bam() -> None:
    # upstream case: picard/addorreplacereadgroups "sarscov2 - bam" (tests/bam.config)
    bam = await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam")
    result = await picard_addorreplacereadgroups(bam, prefix="test.replaced", args="--CREATE_INDEX " + READ_GROUP_ARGS)
    assert result.reads.name == "test.replaced.bam", result.reads.name
    await assert_reads_md5(result.reads, "894549ee3ced6b5ca2eed2563a985217", label="read-grouped bam")
    assert result.bai is not None and result.bai.name == "test.replaced.bai", result.bai
    header = await text_lines(await samtools_view(bam=result.reads, args="-H"))
    read_groups = [set(h.split("\t")[1:]) for h in header if h.startswith("@RG")]
    assert read_groups == [{"ID:1", "LB:test", "PL:ILLUMINA", "PU:bc1", "SM:test"}], read_groups


@env.task
async def test_picard_addorreplacereadgroups_cram() -> None:
    # upstream case: picard/addorreplacereadgroups "homo_sapiens - cram" (tests/cram.config)
    cram = await fixture(HUMAN + "illumina/cram/test.paired_end.sorted.cram")
    fasta = await fixture(HUMAN + "genome/genome.fasta")
    fai = await fixture(HUMAN + "genome/genome.fasta.fai")
    result = await picard_addorreplacereadgroups(
        cram, prefix="test.replaced", fasta=fasta, fai=fai, args=READ_GROUP_ARGS
    )
    assert result.reads.name == "test.replaced.cram", result.reads.name
    assert result.bai is None, result.bai
    # As nft-bam's getReadsMD5 over the decoded CRAM: one md5 of every record's SEQ.
    sam = await text_lines(await samtools_view(bam=result.reads, fasta=fasta, fai=fai))
    actual = hashlib.md5("".join(line.split("\t")[9] for line in sam).encode()).hexdigest().lstrip("0")
    assert actual == "2f11e4fe3390b8ad0a1852616fd1da04", f"reads md5 {actual}"


@env.task
async def test_picard_createsequencedictionary() -> None:
    # upstream case: picard/createsequencedictionary "sarscov2 - fasta"
    fasta = await fixture(SC2 + "genome/genome.fasta")
    out = await picard_createsequencedictionary(fasta, prefix="test")
    assert out.name == "test.dict", out.name
    expected = await text_lines(await fixture(SC2 + "genome/genome.dict"))

    def records(lines: list[str]) -> list[str]:
        return ["\t".join(f for f in line.split("\t") if not f.startswith("UR:")) for line in lines]

    assert records(await text_lines(out)) == records(expected), await text_lines(out)


async def check_collecthsmetrics(metrics: File, size: int, lines_md5: str) -> None:
    assert metrics.name == "test.CollectHsMetrics.coverage_metrics", metrics.name
    async with metrics.open("rb") as fh:
        data = bytes(await fh.read())
    lines = data.decode().splitlines()[:101]
    del lines[3]  # the date, as upstream drops it
    actual = hashlib.md5("\n".join(lines).encode()).hexdigest()
    assert actual == lines_md5, f"first lines md5 {actual}; header {lines[:2]}"
    assert len(data) == size, f"size {len(data)}, expected {size}"


async def hsmetrics_inputs(baits: str, targets: str) -> tuple[File, File, File, File]:
    bam, bai = await sorted_bam()
    return bam, bai, await fixture(SC2 + baits), await fixture(SC2 + targets)


async def sarscov2_reference(gzipped: bool = False) -> dict[str, File | None]:
    fasta = SC2 + ("genome/genome.fasta.gz" if gzipped else "genome/genome.fasta")
    return {
        "fasta": await fixture(fasta),
        "fai": await fixture(fasta + ".fai"),
        "sequence_dict": await fixture(SC2 + "genome/genome.dict"),
        "gzi": await fixture(fasta + ".gzi") if gzipped else None,
    }


INTERVAL_LISTS = ("genome/picard/baits.interval_list", "genome/picard/targets.interval_list")


@env.task
async def test_picard_collecthsmetrics_bam() -> None:
    # upstream case: picard/collecthsmetrics "sarscov2 - bam"
    inputs = await hsmetrics_inputs(*INTERVAL_LISTS)
    metrics = await picard_collecthsmetrics(*inputs, prefix="test", **await sarscov2_reference())
    await check_collecthsmetrics(metrics, 3598, "d5cad9bf2de71bd7799682b4d5747238")


@env.task
async def test_picard_collecthsmetrics_gzippedfa() -> None:
    # upstream case: picard/collecthsmetrics "sarscov2 - bam - gzippedfa"
    inputs = await hsmetrics_inputs(*INTERVAL_LISTS)
    metrics = await picard_collecthsmetrics(*inputs, prefix="test", **await sarscov2_reference(gzipped=True))
    await check_collecthsmetrics(metrics, 3601, "aae9699e1c2ca98e5e6db5319384021f")


@env.task
async def test_picard_collecthsmetrics_nofasta() -> None:
    # upstream case: picard/collecthsmetrics "sarscov2 - bam - nofasta"
    inputs = await hsmetrics_inputs(*INTERVAL_LISTS)
    metrics = await picard_collecthsmetrics(*inputs, prefix="test")
    await check_collecthsmetrics(metrics, 3548, "0f5deac09c6ca74a4d6bfb6d079365e3")


@env.task
async def test_picard_collecthsmetrics_bed() -> None:
    # upstream case: picard/collecthsmetrics "sarscov2 - bam - bed"
    inputs = await hsmetrics_inputs("genome/bed/baits.bed", "genome/bed/test.bed")
    metrics = await picard_collecthsmetrics(*inputs, prefix="test", **await sarscov2_reference())
    await check_collecthsmetrics(metrics, 3595, "c57c8c73ed51b3fc9d7070a6758d06fc")


@env.task
async def test_picard_collecthsmetrics_samebed() -> None:
    # upstream case: picard/collecthsmetrics "sarscov2 - bam - samebed"
    inputs = await hsmetrics_inputs("genome/bed/baits.bed", "genome/bed/baits.bed")
    metrics = await picard_collecthsmetrics(*inputs, prefix="test", **await sarscov2_reference())
    await check_collecthsmetrics(metrics, 3586, "a9162e3220034a5a2aba90e9d7d7e945")


@env.task
async def test_picard_mergesamfiles_bams() -> None:
    # upstream case: picard/mergesamfiles "test-picard-mergesamfiles"
    bams = [
        await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        await fixture(SC2 + "illumina/bam/test.single_end.bam"),
    ]
    merged = await picard_mergesamfiles(bams, prefix="test")
    assert merged.name == "test.bam", merged.name
    await assert_reads_md5(merged, "730da54c088350ff625d34b95e623dca", label="merged bam")


@env.task
async def test_picard_sortsam_queryname() -> None:
    # upstream case: picard/sortsam "test-picard-sortsam"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.single_end.bam")
    sorted_bam = await picard_sortsam(bam, "queryname", prefix="test.sorted")
    assert sorted_bam.name == "test.sorted.bam", f"wrote {sorted_bam.name}"
    await assert_reads_md5(sorted_bam, "b9847fed94d2b7286e18caaa099658ce", label="sorted bam")


@env.task
async def test_picard_filtersamreads() -> None:
    # upstream cases: picard/filtersamreads "sarscov2 - bam" (after SortSam by queryname), "sarscov2 - bam -
    # readlist" and "... - readlist - fasta" (ext.prefix "<id>.filtered", ext.args "--CREATE_INDEX true")
    bam, readlist, fasta = await asyncio.gather(
        fixture(SC2 + "illumina/bam/test.single_end.sorted.bam"),
        fixture(SC2 + "illumina/picard/test.single_end.bam.readlist.txt"),
        fixture(SC2 + "genome/genome.fasta"),
    )
    queryname_sorted = await picard_sortsam(bam, "queryname", prefix="test.sorted")
    common = {"prefix": "test.filtered", "args": "--CREATE_INDEX true"}
    aligned, listed, listed_fasta = await asyncio.gather(
        picard_filtersamreads(queryname_sorted, "includeAligned", **common),
        picard_filtersamreads(bam, "includeReadList", readlist=readlist, **common),
        picard_filtersamreads(bam, "includeReadList", readlist=readlist, fasta=fasta, **common),
    )
    for result, label, md5, indexed in (
        (aligned, "bam", "b9847fed94d2b7286e18caaa099658ce", False),
        (listed, "readlist", "3c2a0fcd376bd70060073cb361fe91c1", True),
        (listed_fasta, "readlist - fasta", "3c2a0fcd376bd70060073cb361fe91c1", True),
    ):
        assert result.bam.name == "test.filtered.bam", f"{label}: {result.bam.name}"
        await assert_reads_md5(result.bam, md5, label=label)
        if indexed:
            assert result.bai is not None and result.bai.name == "test.filtered.bai", f"{label}: {result.bai}"
        else:
            assert result.bai is None, f"{label}: unexpected {result.bai}"


async def check_collectwgsmetrics(metrics: File) -> None:
    assert metrics.name == "test.CollectWgsMetrics.coverage_metrics", metrics.name
    async with metrics.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    # Upstream's whole snapshot: this (space-separated) string is not in the tab-separated table.
    assert ("coverage    high_quality_coverage_count" in text) is False, "unexpected space-separated header"
    lines = text.splitlines()
    assert lines[1].startswith("# CollectWgsMetrics ") and "--REFERENCE_SEQUENCE genome.fasta" in lines[1], lines[1]
    assert "## METRICS CLASS\tpicard.analysis.WgsMetrics" in lines, "no WGS metrics table"
    histogram = lines.index("## HISTOGRAM\tjava.lang.Integer") if "## HISTOGRAM\tjava.lang.Integer" in lines else -1
    assert histogram > 0 and lines[histogram + 1].split("\t")[:2] == ["coverage", "high_quality_coverage_count"], (
        f"no coverage histogram: {lines[histogram : histogram + 2]}"
    )


@env.task
async def test_picard_collectwgsmetrics() -> None:
    # upstream case: picard/collectwgsmetrics "test-picard-collectwgsmetrics"
    bam, bai = await sorted_bam()
    fasta = await fixture(SC2 + "genome/genome.fasta")
    fai = await fixture(SC2 + "genome/genome.fasta.fai")
    await check_collectwgsmetrics(await picard_collectwgsmetrics(bam, fasta, prefix="test", bai=bai, fai=fai))


@env.task
async def test_picard_collectwgsmetrics_with_interval() -> None:
    # upstream case: picard/collectwgsmetrics "test-picard-collectwgsmetrics-with-interval" (no BAM index)
    bam = await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam")
    fasta = await fixture(SC2 + "genome/genome.fasta")
    fai = await fixture(SC2 + "genome/genome.fasta.fai")
    intervals = await fixture(SC2 + "genome/picard/baits.interval_list")
    metrics = await picard_collectwgsmetrics(bam, fasta, prefix="test", fai=fai, intervals=intervals)
    await check_collectwgsmetrics(metrics)
    header = (await text_lines(metrics))[1]
    assert "--INTERVALS baits.interval_list" in header, header


INSERT_SIZE_HEADER = (
    "MEDIAN_INSERT_SIZE\tMODE_INSERT_SIZE\tMEDIAN_ABSOLUTE_DEVIATION\tMIN_INSERT_SIZE\tMAX_INSERT_SIZE\t"
    "MEAN_INSERT_SIZE\tSTANDARD_DEVIATION\tREAD_PAIRS\tPAIR_ORIENTATION\tWIDTH_OF_10_PERCENT\tWIDTH_OF_20_PERCENT\t"
    "WIDTH_OF_30_PERCENT\tWIDTH_OF_40_PERCENT\tWIDTH_OF_50_PERCENT\tWIDTH_OF_60_PERCENT\tWIDTH_OF_70_PERCENT\t"
    "WIDTH_OF_80_PERCENT\tWIDTH_OF_90_PERCENT\tWIDTH_OF_95_PERCENT\tWIDTH_OF_99_PERCENT\tSAMPLE\tLIBRARY\tREAD_GROUP"
)


@env.task
async def test_picard_collectinsertsizemetrics() -> None:
    # upstream case: picard/collectinsertsizemetrics "test-picard-collectinsertsizemetrics"
    bam = await fixture(SC2 + "illumina/bam/test.paired_end.bam")
    result = await picard_collectinsertsizemetrics(bam, prefix="test")
    expected = [
        "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics",
        INSERT_SIZE_HEADER,
        "209\t159\t46\t77\t364\t207.659794\t66.769018\t97\tFR\t25\t49\t59\t77\t93\t123\t145\t183\t223\t255\t311\t\t\t",
        "",
    ]
    lines = (await text_lines(result.metrics))[5:9]
    assert lines == expected, lines
    assert result.histogram.name == "test.pdf", result.histogram.name
    async with result.histogram.open("rb") as fh:
        assert bytes(await fh.read(5)) == b"%PDF-", "histogram is not a PDF"


@env.task
async def test_picard_collectrnaseqmetrics() -> None:
    # upstream case: picard/collectrnaseqmetrics "sarscov2 - fasta - gtf" (refFlat from ucsc/gtftogenepred)
    gtf = await fixture(SC2 + "genome/genome.gtf")
    genepred = await gtftogenepred(gtf, prefix="test", args="-genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons")
    assert genepred.refflat is not None
    bam = await fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam")
    fasta = await fixture(SC2 + "genome/genome.fasta")
    result = await picard_collectrnaseqmetrics(
        bam,
        genepred.refflat,
        prefix="test",
        fasta=fasta,
        args="--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND",
    )
    assert result.metrics.name == "test.rna_metrics", result.metrics.name
    text = "\n".join(await text_lines(result.metrics))
    assert "CollectRnaSeqMetrics" in text, "metrics don't name CollectRnaSeqMetrics"
    assert result.pdfs == [], [f.name for f in result.pdfs]


@env.task
async def test_picard_bedtointervallist() -> None:
    # upstream case: picard/bedtointervallist "test-picard-bedtointervallist"
    bed = await fixture(SC2 + "genome/bed/test.bed")
    sequence_dict = await fixture(SC2 + "genome/genome.dict")
    intervals = await picard_bedtointervallist(bed, sequence_dict, prefix="test")
    assert intervals.name == "test.intervallist", intervals.name
    await assert_md5(intervals, "e51101c9357fb2d59fd30e370eefa39c", label="interval list")


async def picard_vcf_records(vcf: File) -> tuple[int, list[str]]:
    """A gzipped VCF's sample count and record lines."""
    async with vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    header = next(x for x in lines if x.startswith("#CHROM"))
    return len(header.split("\t")) - 9, [x for x in lines if not x.startswith("#")]


@env.task
async def test_picard_liftovervcf() -> None:
    hs = "genomics/homo_sapiens/"
    vcf, sequence_dict, fasta, chain = await asyncio.gather(
        fixture(hs + "illumina/gvcf/test.genome.vcf"),
        fixture(hs + "genome/genome.dict"),
        fixture(hs + "genome/genome.fasta"),
        fixture(hs + "genome/genome.chain.gz"),
    )
    result = await picard_liftovervcf(
        vcf, sequence_dict, fasta, chain, prefix="test", args="--WARN_ON_MISSING_CONTIG true"
    )
    assert result.lifted.name == "test.lifted.vcf.gz", result.lifted.name
    samples, records = await picard_vcf_records(result.lifted)
    assert (samples, len(records)) == (1, 0), f"lifted VCF: {samples} samples, {len(records)} records"
    samples, records = await picard_vcf_records(result.unlifted)
    md5 = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
    assert md5 == "39a9de5185d94289283bd27cfcdeba97", f"unlifted variantsMD5 {md5}"


tests = [
    test_markduplicates_sorted_bam,
    test_picard_collectmultiplemetrics,
    test_picard_collectmultiplemetrics_nofasta,
    test_picard_collectmultiplemetrics_cram,
    test_picard_addorreplacereadgroups_bam,
    test_picard_addorreplacereadgroups_cram,
    test_picard_createsequencedictionary,
    test_picard_collecthsmetrics_bam,
    test_picard_collecthsmetrics_gzippedfa,
    test_picard_collecthsmetrics_nofasta,
    test_picard_collecthsmetrics_bed,
    test_picard_collecthsmetrics_samebed,
    test_picard_mergesamfiles_bams,
    test_picard_sortsam_queryname,
    test_picard_filtersamreads,
    test_picard_collectwgsmetrics,
    test_picard_collectwgsmetrics_with_interval,
    test_picard_collectinsertsizemetrics,
    test_picard_collectrnaseqmetrics,
    test_picard_bedtointervallist,
    test_picard_liftovervcf,
]
