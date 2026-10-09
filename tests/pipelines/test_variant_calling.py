"""End-to-end tests for flyte_bio.pipelines.variant_calling on upstream's test genome.

Upstream's test profile (``testdata.nf-core.sarek``) uses a 40 kb slice of
chr22 with prebuilt reference files, so references built from the FASTA can
be checked against them.
"""

import asyncio
import gzip
import hashlib

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.samtools import samtools_view
from flyte_bio.pipelines.variant_calling import (
    SampleQC,
    call_variants,
    map_reads,
    mark_duplicates,
    multiqc_report,
    prepare_reference,
    recalibrate,
    samples_from_samplesheet,
    variant_calling,
    vcf_qc,
)
from flyte_bio.samplesheet import read_samplesheet
from tests.framework import assert_md5, assert_nonempty, env, fixture, stage_samplesheet

GENOME = "genomics/homo_sapiens/genome/"
FASTA = GENOME + "genome.fasta"
FAI = GENOME + "genome.fasta.fai"
DICT = GENOME + "genome.dict"
INTERVALS = GENOME + "genome.interval_list"
DBSNP = GENOME + "vcf/dbsnp_146.hg38.vcf.gz"
KNOWN_INDELS = GENOME + "vcf/mills_and_1000G.indels.vcf.gz"
MULTI_INTERVALS = GENOME + "genome.multi_intervals.bed"  # two intervals: two scatter chunks at 20 nt/s
MAPPED_BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"  # upstream's mapped_single_bam.csv


async def vcf_records(vcf: File) -> tuple[list[str], list[str]]:
    """A gzipped VCF's (header lines, record lines)."""
    async with vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    return [x for x in lines if x.startswith("#")], [x for x in lines if not x.startswith("#")]


CSV = "https://raw.githubusercontent.com/nf-core/sarek/3.10.0/tests/csv/3.0/"
FASTQ_SAMPLESHEET = CSV + "fastq_single.csv"  # one sample, two lanes (the same FASTQ pair)


async def text(f: File) -> str:
    async with f.open("rb") as fh:
        return bytes(await fh.read()).decode()


def sequences(dictionary: str) -> list[str]:
    """A sequence dictionary's @SQ lines without the UR field (where the FASTA was)."""
    return [
        "\t".join(field for field in line.split("\t") if not field.startswith("UR:"))
        for line in dictionary.splitlines()
        if line.startswith("@SQ")
    ]


@env.task
async def test_prepare_reference_from_fasta() -> None:
    fasta, upstream_fai, upstream_dict = await asyncio.gather(fixture(FASTA), fixture(FAI), fixture(DICT))
    ref = await prepare_reference(fasta)

    # The FASTA index is a pure function of the FASTA: upstream's prebuilt one matches byte for byte.
    async with upstream_fai.open("rb") as fh:
        expected_fai = hashlib.md5(bytes(await fh.read())).hexdigest()
    await assert_md5(ref.fai, expected_fai, label="built .fai")
    assert sequences(await text(ref.sequence_dict)) == sequences(await text(upstream_dict))
    assert ref.sequence_dict.path.endswith("/genome.dict"), ref.sequence_dict.path
    for ext in ("amb", "ann", "bwt", "pac", "sa"):
        f = await ref.bwa_index.get_file(f"genome.{ext}")
        assert f is not None, f"bwa index has no genome.{ext}"
    # No intervals given: one per sequence, named after the FASTA.
    assert ref.intervals is not None and ref.intervals.path.endswith("/genome.bed"), ref.intervals
    assert await text(ref.intervals) == "chr22\t0\t40001\n", await text(ref.intervals)
    assert [f.path.rsplit("/", 1)[-1] for f in ref.intervals_split] == ["chr22_1-40001.bed"]
    assert ref.dbsnp is None and ref.known_indels == []


@env.task
async def test_prepare_reference_test_profile() -> None:
    # As upstream's test profile supplies it: fai, dict and an interval list, known sites without indexes.
    fasta, fai, dict_file, intervals, dbsnp, indels = await asyncio.gather(
        fixture(FASTA), fixture(FAI), fixture(DICT), fixture(INTERVALS), fixture(DBSNP), fixture(KNOWN_INDELS)
    )
    ref = await prepare_reference(
        fasta,
        fai=fai,
        sequence_dict=dict_file,
        dbsnp=dbsnp,
        known_indels=[indels],
        intervals=intervals,
        nucleotides_per_second=20,
    )
    assert ref.fai is fai and ref.sequence_dict is dict_file
    # The interval list converted to BED (GATK may add name/strand columns).
    assert ref.intervals is not None and ref.intervals.path.endswith("/genome.bed"), ref.intervals
    rows = [line.split("\t")[:3] for line in (await text(ref.intervals)).splitlines()]
    assert rows == [["chr22", "0", "40001"]], rows
    assert [f.path.rsplit("/", 1)[-1] for f in ref.intervals_split] == ["chr22_1-40001.bed"]
    assert ref.dbsnp is not None and ref.dbsnp.path.endswith("/dbsnp_146.hg38.vcf.gz"), ref.dbsnp
    assert ref.dbsnp_tbi is not None and ref.dbsnp_tbi.path.endswith("/dbsnp_146.hg38.vcf.gz.tbi")
    assert [f.path.rsplit("/", 1)[-1] for f in ref.known_indels] == ["mills_and_1000G.indels.vcf.gz"]
    assert [f.path.rsplit("/", 1)[-1] for f in ref.known_indels_tbi] == ["mills_and_1000G.indels.vcf.gz.tbi"]
    for f in (ref.dbsnp_tbi, *ref.known_indels_tbi):
        await assert_nonempty(f, label=f.path.rsplit("/", 1)[-1])


@env.task
async def test_map_reads() -> None:
    sheet = await stage_samplesheet(FASTQ_SAMPLESHEET, ("fastq_1", "fastq_2"))
    (sample,) = samples_from_samplesheet(await read_samplesheet(sheet))
    assert (sample.patient, sample.sex, sample.status, sample.id) == ("test", "XX", 0, "test"), sample
    assert [lane.id for lane in sample.lanes] == ["test_L1", "test_L2"]
    assert all(lane.fastq_2 is not None for lane in sample.lanes)

    fasta, fai, dict_file = await asyncio.gather(fixture(FASTA), fixture(FAI), fixture(DICT))
    ref = await prepare_reference(fasta, fai=fai, sequence_dict=dict_file)
    mapped = await map_reads(sample, ref)

    assert [b.path.rsplit("/", 1)[-1] for b in mapped.lane_bams] == [
        "test-test_L1.sorted.bam",
        "test-test_L2.sorted.bam",
    ]
    # Upstream's read groups: the test reads' names carry no flowcell, so ID is <sample>.<lane>.
    for lane, rg in zip(("test_L1", "test_L2"), mapped.read_groups):
        assert rg == (f"@RG\\tID:test.{lane}\\tPU:{lane}\\tSM:test_test\\tLB:test\\tDS:{fasta.path}\\tPL:ILLUMINA"), rg
    for lane, bam in zip(("test_L1", "test_L2"), mapped.lane_bams):
        header = await text(await samtools_view(bam=bam, args="-H"))
        rg_lines = [line for line in header.splitlines() if line.startswith("@RG")]
        assert rg_lines and rg_lines[0].startswith(f"@RG\tID:test.{lane}\tPU:{lane}\tSM:test_test"), rg_lines
        assert "SO:coordinate" in header.splitlines()[0], header.splitlines()[0]
        await assert_nonempty(bam, label=f"{lane} bam")
    # Both lanes are the same FASTQ pair, so they hold the same alignments.
    counts = [len((await text(await samtools_view(bam=b))).splitlines()) for b in mapped.lane_bams]
    assert counts[0] == counts[1] > 0, counts
    assert len(mapped.fastqc) == 2
    for report_dir, lane in zip(mapped.fastqc, ("test_L1", "test_L2")):
        assert await report_dir.get_file(f"test-{lane}_1_fastqc.html") is not None, f"no FastQC report for {lane}"


def duplication_rate(metrics: str) -> float:
    """PERCENT_DUPLICATION of a Picard duplication-metrics file (single library)."""
    lines = metrics.splitlines()
    i = next(n for n, line in enumerate(lines) if line.startswith("## METRICS CLASS"))
    header, values = lines[i + 1].split("\t"), lines[i + 2].split("\t")
    return float(values[header.index("PERCENT_DUPLICATION")])


@env.task
async def test_mark_duplicates() -> None:
    sheet = await stage_samplesheet(FASTQ_SAMPLESHEET, ("fastq_1", "fastq_2"))
    (sample,) = samples_from_samplesheet(await read_samplesheet(sheet))
    fasta, fai, dict_file = await asyncio.gather(fixture(FASTA), fixture(FAI), fixture(DICT))
    ref = await prepare_reference(fasta, fai=fai, sequence_dict=dict_file)
    marked = await mark_duplicates(await map_reads(sample, ref), ref)

    assert marked.cram.path.endswith("/test.md.cram"), marked.cram.path
    assert marked.crai.path.endswith("/test.md.cram.crai"), marked.crai.path
    assert marked.metrics.path.endswith("/test.md.cram.metrics"), marked.metrics.path
    # The two lanes are the same FASTQ pair, so at least half the reads are duplicates.
    rate = duplication_rate(await text(marked.metrics))
    assert rate >= 0.5, rate
    # One CRAM holding both lanes' read groups.
    header = await text(await samtools_view(bam=marked.cram, args="-H"))
    ids = sorted(line.split("\t")[1] for line in header.splitlines() if line.startswith("@RG"))
    assert ids == ["ID:test.test_L1", "ID:test.test_L2"], ids
    stats = await text(marked.qc.samtools_stats)
    assert stats.startswith("# This file was produced by samtools stats"), stats[:200]
    for name in ("test.md.mosdepth.global.dist.txt", "test.md.mosdepth.summary.txt", "test.md.regions.bed.gz"):
        f = await marked.qc.mosdepth.get_file(name)
        assert f is not None, f"mosdepth wrote no {name}"
        await assert_nonempty(f, label=name)


@env.task
async def test_recalibrate() -> None:
    sheet = await stage_samplesheet(FASTQ_SAMPLESHEET, ("fastq_1", "fastq_2"))
    (sample,) = samples_from_samplesheet(await read_samplesheet(sheet))
    fasta, fai, dict_file, dbsnp, indels, intervals = await asyncio.gather(
        fixture(FASTA), fixture(FAI), fixture(DICT), fixture(DBSNP), fixture(KNOWN_INDELS), fixture(MULTI_INTERVALS)
    )
    ref = await prepare_reference(
        fasta,
        fai=fai,
        sequence_dict=dict_file,
        dbsnp=dbsnp,
        known_indels=[indels],
        intervals=intervals,
        nucleotides_per_second=20,
    )
    assert len(ref.intervals_split) == 2, ref.intervals_split
    marked = await mark_duplicates(await map_reads(sample, ref), ref)
    recal = await recalibrate(marked, sample.id, ref)

    # Scattered over two chunks: the tables are gathered and the CRAMs merged under upstream's names.
    assert recal.table.path.endswith("/test.recal.table"), recal.table.path
    assert recal.cram.path.endswith("/test.recal.cram"), recal.cram.path
    table = await text(recal.table)
    assert table.startswith("#:GATKReport"), table[:100]
    assert "RecalTable" in table, "no recalibration tables in the report"
    header = await text(await samtools_view(bam=recal.cram, args="-H"))
    assert any(line.startswith("@PG") and "ApplyBQSR" in line for line in header.splitlines()), "no ApplyBQSR @PG"
    records = (await text(await samtools_view(bam=recal.cram, args="-c"))).strip()
    assert int(records) > 0, records
    await assert_nonempty(recal.crai, label="test.recal.cram.crai")
    for name in ("test.recal.mosdepth.global.dist.txt", "test.recal.mosdepth.summary.txt"):
        assert await recal.qc.mosdepth.get_file(name) is not None, f"mosdepth wrote no {name}"
    assert (await text(recal.qc.samtools_stats)).startswith("# This file was produced by samtools stats")


@env.task
async def test_call_variants() -> None:
    # upstream case: tests/variant_calling_haplotypecaller.nf.test, first scenario
    # (--step variant_calling from mapped_single_bam.csv, test profile, --nucleotides_per_second 20)
    fasta, fai, dict_file, dbsnp, indels, intervals, bam, bai = await asyncio.gather(
        fixture(FASTA),
        fixture(FAI),
        fixture(DICT),
        fixture(DBSNP),
        fixture(KNOWN_INDELS),
        fixture(INTERVALS),
        fixture(MAPPED_BAM),
        fixture(MAPPED_BAM + ".bai"),
    )
    ref = await prepare_reference(
        fasta,
        fai=fai,
        sequence_dict=dict_file,
        dbsnp=dbsnp,
        known_indels=[indels],
        intervals=intervals,
        nucleotides_per_second=20,
    )
    # The test profile's filter arguments.
    result = await call_variants(bam, bai, "test", ref, filter_args="--info-key CNN_1D --indel-tranche 0")
    assert result.calls.vcf.path.endswith("/test.haplotypecaller.vcf.gz"), result.calls.vcf.path
    assert result.scored is not None and result.scored.vcf.path.endswith("/test.cnn.vcf.gz"), result.scored
    assert result.filtered is not None, result
    assert result.filtered.vcf.path.endswith("/test.haplotypecaller.filtered.vcf.gz"), result.filtered.vcf.path
    assert result.final is result.filtered
    header, calls = await vcf_records(result.calls.vcf)
    assert calls, "no calls"
    assert any(line.startswith("##GATKCommandLine=<ID=HaplotypeCaller") and "CONSERVATIVE" in line for line in header)
    # dbSNP IDs are annotated where the calls hit known sites.
    filtered_header, filtered = await vcf_records(result.filtered.vcf)
    assert len(filtered) == len(calls), (len(filtered), len(calls))
    assert any(h.startswith("##FILTER=<ID=CNN_1D_") for h in filtered_header), "no tranche filters"
    assert all("CNN_1D=" in r.split("\t")[7] for r in filtered), "records without a CNN score"

    # Scattered over two chunks (and unfiltered): the per-chunk calls are merged in genome order.
    multi = await fixture(MULTI_INTERVALS)
    scattered_ref = await prepare_reference(
        fasta, fai=fai, sequence_dict=dict_file, dbsnp=dbsnp, intervals=multi, nucleotides_per_second=20
    )
    scattered = await call_variants(bam, bai, "test", scattered_ref, skip_filter=True)
    assert scattered.filtered is None and scattered.final is scattered.calls
    assert scattered.calls.vcf.path.endswith("/test.haplotypecaller.vcf.gz"), scattered.calls.vcf.path
    _, merged = await vcf_records(scattered.calls.vcf)
    positions = [int(r.split("\t")[1]) for r in merged]
    assert merged and positions == sorted(positions), positions[:10]


@env.task
async def test_fastq_to_report() -> None:
    # FASTQ -> MarkDuplicates -> BQSR -> HaplotypeCaller -> VCF QC -> MultiQC, on the test profile's
    # reference (one interval chunk at 20 nt/s).
    sheet = await stage_samplesheet(FASTQ_SAMPLESHEET, ("fastq_1", "fastq_2"))
    (sample,) = samples_from_samplesheet(await read_samplesheet(sheet))
    fasta, fai, dict_file, dbsnp, indels, intervals = await asyncio.gather(
        fixture(FASTA), fixture(FAI), fixture(DICT), fixture(DBSNP), fixture(KNOWN_INDELS), fixture(INTERVALS)
    )
    ref = await prepare_reference(
        fasta,
        fai=fai,
        sequence_dict=dict_file,
        dbsnp=dbsnp,
        known_indels=[indels],
        intervals=intervals,
        nucleotides_per_second=20,
    )
    mapped = await map_reads(sample, ref)
    marked = await mark_duplicates(mapped, ref)
    recal = await recalibrate(marked, sample.id, ref)
    calls = await call_variants(
        recal.cram, recal.crai, sample.id, ref, filter_args="--info-key CNN_1D --indel-tranche 0"
    )
    # The test reads come from the whole genome but the reference is 40 kb of chr22: the reads that map
    # carry ~8% mismatches, BQSR rightly lowers their qualities, and HaplotypeCaller calls nothing. With no
    # records to score, the CNN filter is skipped (GATK would refuse the VCF) and QC runs on the calls.
    _, records = await vcf_records(calls.calls.vcf)
    assert records == [] and calls.filtered is None and calls.final is calls.calls, (len(records), calls)
    qc = await vcf_qc(calls.final)
    assert qc.bcftools_stats.path.endswith("/test.haplotypecaller.bcftools_stats.txt"), qc.bcftools_stats
    for f, suffix in (
        (qc.tstv_count, "TsTv.count"),
        (qc.tstv_qual, "TsTv.qual"),
        (qc.filter_summary, "FILTER.summary"),
    ):
        assert f.path.endswith(f"/test.haplotypecaller.{suffix}"), f.path
        await assert_nonempty(f, label=suffix)

    report = await multiqc_report(
        [SampleQC(sample.id, fastqc=mapped.fastqc, markduplicates=marked, recalibrated=recal, vcf=qc)]
    )
    await assert_nonempty(report.report, label="multiqc_report.html")
    data = [f.path.rsplit("/", 1)[-1].lower() async for f in report.results.walk() if "_data/" in f.path]
    for module in ("fastqc", "picard", "samtools", "mosdepth", "gatk", "bcftools", "vcftools"):
        assert any(module in n for n in data), f"MultiQC parsed no {module} outputs: {sorted(data)}"


async def published(tree: Dir) -> set[str]:
    base = tree.path.rstrip("/")
    return {f.path[len(base) :].lstrip("/") async for f in tree.walk()}


async def profile_inputs() -> dict:
    """variant_calling() keyword arguments for upstream's test profile (one interval chunk at 20 nt/s)."""
    sheet, fasta, fai, dict_file, dbsnp, indels, intervals = await asyncio.gather(
        stage_samplesheet(FASTQ_SAMPLESHEET, ("fastq_1", "fastq_2")),
        fixture(FASTA),
        fixture(FAI),
        fixture(DICT),
        fixture(DBSNP),
        fixture(KNOWN_INDELS),
        fixture(INTERVALS),
    )
    return {
        "samplesheet": sheet,
        "fasta": fasta,
        "fai": fai,
        "sequence_dict": dict_file,
        "dbsnp": dbsnp,
        "known_indels": [indels],
        "intervals": intervals,
        "nucleotides_per_second": 20,
        "haplotypecaller_filter_args": "--info-key CNN_1D --indel-tranche 0",
    }


@env.task
async def test_variant_calling() -> None:
    # The whole germline route from upstream's FASTQ samplesheet, published in upstream's --outdir layout.
    result = await variant_calling(**await profile_inputs(), publish_results=True)
    (r,) = result.samples
    assert r.sample.id == "test" and r.recalibrated is not None and r.vcf_qc is not None
    assert result.multiqc is not None and result.outdir is not None
    paths = await published(result.outdir)
    for path in [
        "reports/fastqc/test-test_L1/test-test_L1_1_fastqc.html",
        "reports/fastqc/test-test_L2/test-test_L2_2_fastqc.html",
        "preprocessing/markduplicates/test/test.md.cram",
        "preprocessing/markduplicates/test/test.md.cram.crai",
        "reports/markduplicates/test/test.md.cram.metrics",
        "reports/samtools/test/test.md.cram.stats",
        "reports/mosdepth/test/test.md.mosdepth.summary.txt",
        "preprocessing/recal_table/test/test.recal.table",
        "preprocessing/recalibrated/test/test.recal.cram",
        "preprocessing/recalibrated/test/test.recal.cram.crai",
        "reports/samtools/test/test.recal.cram.stats",
        "reports/mosdepth/test/test.recal.mosdepth.summary.txt",
        "variant_calling/haplotypecaller/test/test.haplotypecaller.vcf.gz",
        "variant_calling/haplotypecaller/test/test.haplotypecaller.vcf.gz.tbi",
        "reports/bcftools/haplotypecaller/test/test.haplotypecaller.bcftools_stats.txt",
        "reports/vcftools/haplotypecaller/test/test.haplotypecaller.TsTv.count",
        "multiqc/multiqc_report.html",
    ]:
        assert path in paths, f"not published: {path}"


@env.task
async def test_variant_calling_skips() -> None:
    # Without duplicate marking and recalibration the lanes are merged into <sample>.sorted.cram and called
    # as they are. Published to an explicit outdir. The filter is off too: these calls include indels that
    # no known-indel site overlaps, which GATK's tranche filter refuses (upstream would fail the same way);
    # test_call_variants covers the filter.
    outdir = flyte.ctx().raw_data_path.get_random_remote_path("variant_calling_outdir")
    result = await variant_calling(
        **await profile_inputs(),
        skip_markduplicates=True,
        skip_baserecalibrator=True,
        skip_haplotypecaller_filter=True,
        outdir=outdir,
    )
    (r,) = result.samples
    assert r.markduplicates.metrics is None and r.recalibrated is None
    assert r.markduplicates.cram.path.endswith("/test.sorted.cram"), r.markduplicates.cram.path
    assert r.calls.filtered is None and r.calls.final is r.calls.calls
    assert result.outdir is not None and result.outdir.path.rstrip("/") == outdir.rstrip("/")
    paths = await published(result.outdir)
    for path in [
        "preprocessing/mapped/test/test.sorted.cram",
        "reports/samtools/test/test.sorted.cram.stats",
        "reports/mosdepth/test/test.sorted.mosdepth.summary.txt",
        "variant_calling/haplotypecaller/test/test.haplotypecaller.vcf.gz",
        "reports/bcftools/haplotypecaller/test/test.haplotypecaller.bcftools_stats.txt",
        "multiqc/multiqc_report.html",
    ]:
        assert path in paths, f"not published: {path}"
    assert not any(p.startswith(("preprocessing/recalibrated", "reports/markduplicates")) for p in paths), paths


tests = [
    test_prepare_reference_from_fasta,
    test_prepare_reference_test_profile,
    test_map_reads,
    test_mark_duplicates,
    test_recalibrate,
    test_call_variants,
    test_fastq_to_report,
    test_variant_calling,
    test_variant_calling_skips,
]
