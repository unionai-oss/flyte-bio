"""variant_calling — germline short-variant calling from FASTQ (BWA-MEM + GATK).

Ports the germline GATK HaplotypeCaller route of the upstream variant-calling
pipeline (3.10.0):

- :func:`prepare_reference` — index the genome (FASTA index, sequence
  dictionary, BWA index), tabix-index the known-sites VCFs, and build the
  calling intervals and their scatter-gather chunks — each only when not
  supplied, as upstream.
- :func:`samples_from_samplesheet` — group samplesheet rows (``patient, sex,
  status, sample, lane, fastq_1, fastq_2``; one per lane) into samples.
- :func:`map_reads` — FastQC and BWA-MEM each lane of a sample, with
  upstream's read groups, into coordinate-sorted BAMs.
- :func:`mark_duplicates` — GATK MarkDuplicates across a sample's lanes into
  one CRAM, then samtools stats and mosdepth on it.
- :func:`recalibrate` — base-quality score recalibration: BaseRecalibrator
  and ApplyBQSR scattered over the interval chunks, gathered into one CRAM,
  then QC'd like the duplicate-marked one.
- :func:`call_variants` — GATK HaplotypeCaller scattered over the interval
  chunks and merged, then (by default) CNNScoreVariants +
  FilterVariantTranches.
- :func:`vcf_qc` — bcftools stats and vcftools Ts/Tv and FILTER summaries of
  the calls.
- :func:`multiqc_report` — one MultiQC report over every sample's QC, with
  upstream's config.
- :func:`variant_calling` — all of the above for the samples in a
  samplesheet, in parallel, optionally publishing the results in upstream's
  ``--outdir`` layout (:func:`variant_calling_layout`).

Upstream splits FASTQs into 50M-read chunks with fastp by default to map
them in parallel; that isn't ported (each lane maps as one job), which
changes how the work is spread, not the alignments.

These are plain async functions, not tasks: they run inside the caller's
task and fan out to the module tasks, so the caller's
:class:`flyte.TaskEnvironment` must ``depends_on`` :data:`flyte_bio.modules.env`.
"""

import asyncio
import tempfile
import urllib.request
import zlib
from dataclasses import dataclass, field
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.bcftools import bcftools_stats
from flyte_bio.modules.bwa import bwa_index as build_bwa_index
from flyte_bio.modules.bwa import bwa_mem
from flyte_bio.modules.fastqc import fastqc
from flyte_bio.modules.gatk4 import (
    IndexedVcf,
    gatk4_applybqsr,
    gatk4_baserecalibrator,
    gatk4_cnnscorevariants,
    gatk4_createsequencedictionary,
    gatk4_filtervarianttranches,
    gatk4_gatherbqsrreports,
    gatk4_haplotypecaller,
    gatk4_intervallisttobed,
    gatk4_markduplicates,
    gatk4_mergevcfs,
)
from flyte_bio.modules.htslib import htslib_bgziptabix
from flyte_bio.modules.intervals import build_intervals, create_intervals_bed
from flyte_bio.modules.mosdepth import mosdepth
from flyte_bio.modules.multiqc import MultiQCResult, multiqc
from flyte_bio.modules.samtools import samtools_faidx, samtools_index, samtools_merge, samtools_stats
from flyte_bio.modules.vcftools import vcftools
from flyte_bio.publish import expand, name, publish
from flyte_bio.samplesheet import Row, read_samplesheet
from flyte_bio.scripts import path as asset


@dataclass
class Reference:
    fasta: File
    fai: File
    sequence_dict: File
    bwa_index: Dir
    dbsnp: File | None = None
    dbsnp_tbi: File | None = None
    known_indels: list[File] = field(default_factory=list)
    known_indels_tbi: list[File] = field(default_factory=list)
    # All calling intervals in one BED; None with no_intervals (whole genome, unsplit).
    intervals: File | None = None
    # The intervals split for scatter-gather, longest-running first; empty with no_intervals.
    intervals_split: list[File] = field(default_factory=list)


def base_name(f: File) -> str:
    """A file's name without its last extension (Nextflow's ``baseName``)."""
    n = name(f)
    return n.rsplit(".", 1)[0] if "." in n else n


async def tabix_index(vcf: File) -> tuple[File, File]:
    """bgzip + tabix a VCF as upstream does: ``<name without last extension>.gz`` and its ``.tbi``."""
    gz, tbi = await htslib_bgziptabix(vcf, prefix=base_name(vcf))
    assert tbi is not None
    return gz, tbi


async def prepare_reference(
    fasta: File,
    fai: File | None = None,
    sequence_dict: File | None = None,
    bwa_index: Dir | None = None,
    dbsnp: File | None = None,
    dbsnp_tbi: File | None = None,
    known_indels: list[File] | None = None,
    known_indels_tbi: list[File] | None = None,
    intervals: File | None = None,
    no_intervals: bool = False,
    nucleotides_per_second: int = 200000,
) -> Reference:
    """Everything the germline route needs from the reference, building what isn't supplied.

    ``intervals`` is a BED or Picard ``.interval_list`` of the regions to call
    (e.g. exome targets); without it the whole genome is used, one interval
    per sequence. The intervals are then split into scatter-gather chunks
    sized by ``nucleotides_per_second`` (or a BED's fifth column of runtime
    estimates). ``no_intervals`` skips intervals altogether: every step runs
    on the whole genome in one go.

    A known-sites VCF supplied without its ``.tbi`` is bgzipped and indexed.
    """
    known_indels = known_indels or []
    if known_indels_tbi and len(known_indels_tbi) != len(known_indels):
        raise ValueError("known_indels_tbi must have one index per known_indels VCF")

    async def get_fai() -> File:
        return fai if fai is not None else (await samtools_faidx(fasta=fasta))[0]

    async def get_dict() -> File:
        return sequence_dict if sequence_dict is not None else await gatk4_createsequencedictionary(fasta)

    async def get_bwa_index() -> Dir:
        return bwa_index if bwa_index is not None else await build_bwa_index(fasta)

    async def get_dbsnp() -> tuple[File | None, File | None]:
        if dbsnp is None or dbsnp_tbi is not None:
            return dbsnp, dbsnp_tbi
        return await tabix_index(dbsnp)

    async def get_known_indels() -> tuple[list[File], list[File]]:
        if known_indels_tbi:
            return known_indels, known_indels_tbi
        indexed = await asyncio.gather(*(tabix_index(vcf) for vcf in known_indels))
        return [gz for gz, _ in indexed], [tbi for _, tbi in indexed]

    ref_fai, ref_dict, index, (dbsnp_vcf, dbsnp_index), (indels, indels_tbi) = await asyncio.gather(
        get_fai(), get_dict(), get_bwa_index(), get_dbsnp(), get_known_indels()
    )

    combined, split = None, []
    if not no_intervals:
        if intervals is None:
            # As upstream: one interval per sequence, named after the FASTA.
            combined = await build_intervals(ref_fai, prefix=base_name(fasta))
            split = await create_intervals_bed(combined, nucleotides_per_second)
        else:
            # As upstream: split the intervals as given; an interval_list is also converted to BED.
            combined = intervals
            if name(intervals).endswith(".interval_list"):
                combined = await gatk4_intervallisttobed(intervals, prefix=base_name(intervals))
            split = await create_intervals_bed(intervals, nucleotides_per_second)

    return Reference(
        fasta=fasta,
        fai=ref_fai,
        sequence_dict=ref_dict,
        bwa_index=index,
        dbsnp=dbsnp_vcf,
        dbsnp_tbi=dbsnp_index,
        known_indels=indels,
        known_indels_tbi=indels_tbi,
        intervals=combined,
        intervals_split=split,
    )


@dataclass
class Lane:
    id: str  # the samplesheet's lane column
    fastq_1: File
    fastq_2: File | None = None  # None for single-end


@dataclass
class Sample:
    """One sample of one patient: its sequencing lanes (mapped separately, merged at duplicate marking)."""

    patient: str
    sex: str
    status: int  # 0 normal, 1 tumor
    id: str  # the samplesheet's sample column
    lanes: list[Lane]


def samples_from_samplesheet(rows: list[Row]) -> list[Sample]:
    """Group samplesheet rows (one per lane) into Samples, in file order, as upstream.

    Upstream's columns for FASTQ input: ``patient``, ``sex``, ``status`` (0
    normal, 1 tumor), ``sample``, ``lane``, ``fastq_1`` and ``fastq_2``
    (empty for single-end). A sample's rows must agree on patient, sex and
    status, and its lanes must be distinct. Only germline (status 0) samples
    are supported so far.
    """
    groups: dict[tuple[str, str], list[Row]] = {}
    for row in rows:
        groups.setdefault((row["patient"], row["sample"]), []).append(row)
    samples = []
    for (patient, sample_id), sample_rows in groups.items():
        lines = [r.line for r in sample_rows]
        for column in ("sex", "status"):
            values = {r[column] for r in sample_rows}
            if len(values) > 1:
                raise ValueError(f"sample {sample_id!r} (lines {lines}): rows disagree on {column} {sorted(values)}")
        status = int(sample_rows[0]["status"] or 0)
        if status != 0:
            raise ValueError(f"sample {sample_id!r} (lines {lines}): tumor samples (status 1) aren't supported yet")
        lanes = [r["lane"] for r in sample_rows]
        if len(lanes) != len(set(lanes)):
            raise ValueError(f"sample {sample_id!r} (lines {lines}): lanes must be distinct, got {lanes}")
        samples.append(
            Sample(
                patient=patient,
                sex=sample_rows[0]["sex"],
                status=status,
                id=sample_id,
                lanes=[
                    Lane(
                        id=r["lane"],
                        fastq_1=File.from_existing_remote(r["fastq_1"]),
                        fastq_2=File.from_existing_remote(r["fastq_2"]) if r["fastq_2"] else None,
                    )
                    for r in sample_rows
                ],
            )
        )
    return samples


async def first_line(fastq: File) -> str:
    """The first line of a gzipped FASTQ, reading only as much as it takes.

    ``https://`` FASTQs are read with the standard library: Flyte's file
    reader can't fetch them reliably inside a task.
    """
    decompressor = zlib.decompressobj(zlib.MAX_WBITS | 16)
    text = b""
    if fastq.path.startswith(("http://", "https://")):

        def read_http() -> bytes:
            data = b""
            with urllib.request.urlopen(fastq.path, timeout=60) as response:
                while b"\n" not in data:
                    chunk = response.read(64 * 1024)
                    if not chunk:
                        break
                    data += decompressor.decompress(chunk)
            return data

        text = await asyncio.to_thread(read_http)
    else:
        async with fastq.open("rb") as fh:
            while b"\n" not in text:
                chunk = bytes(await fh.read(64 * 1024))
                if not chunk:
                    break
                text += decompressor.decompress(chunk)
    return text.split(b"\n", 1)[0].decode(errors="replace")


def flowcell_from_header(header: str) -> str | None:
    """Upstream's flowcell ID from an Illumina read name, or None if it isn't one.

    ``@<instrument>:<lane>:<tile>:<x>:<y>`` gives the instrument;
    ``@<instrument>:<run>:<flowcell>:<lane>:...`` (CASAVA 1.8+) the flowcell.
    """
    fields = header.split(":") if header else []
    if len(fields) == 5:
        return fields[0][1:]
    if len(fields) >= 7:
        return fields[2]
    return None


async def read_group(
    sample: Sample, lane: Lane, fasta: File, seq_platform: str = "ILLUMINA", seq_center: str = ""
) -> tuple[str, str]:
    """Upstream's read group for a lane, and its ``sample_lane_id``: ``(@RG line, id)``.

    The ID is ``<flowcell>.<sample>.<lane>`` when the FASTQ's read names give a
    flowcell (both mates must agree), else ``<sample>.<lane>``.
    """
    flowcell = flowcell_from_header(await first_line(lane.fastq_1))
    if flowcell and lane.fastq_2 is not None and flowcell != flowcell_from_header(await first_line(lane.fastq_2)):
        raise ValueError(f"sample {sample.id!r} lane {lane.id!r}: flowcell IDs of the two mates differ")
    sample_lane_id = f"{flowcell}.{sample.id}.{lane.id}" if flowcell else f"{sample.id}.{lane.id}"
    center = f"CN:{seq_center}\\t" if seq_center else ""
    rg = (
        f"@RG\\tID:{sample_lane_id}\\t{center}PU:{lane.id}\\tSM:{sample.patient}_{sample.id}"
        f"\\tLB:{sample.id}\\tDS:{fasta.path}\\tPL:{seq_platform}"
    )
    return rg, sample_lane_id


@dataclass
class MappedSample:
    sample: Sample
    lane_bams: list[File]  # one coordinate-sorted BAM per lane: <sample>-<lane>.sorted.bam
    read_groups: list[str]
    fastqc: list[Dir]  # FastQC reports per lane (empty when skipped)


async def map_reads(
    sample: Sample,
    reference: Reference,
    seq_platform: str = "ILLUMINA",
    seq_center: str = "",
    skip_fastqc: bool = False,
) -> MappedSample:
    """FastQC and map each lane of ``sample`` with BWA-MEM, as upstream.

    BWA-MEM runs with upstream's ``-K 100000000 -Y`` and the lane's read
    group, piped into ``samtools sort``. Lanes stay separate: MarkDuplicates
    takes them all at once.
    """
    for lane in sample.lanes:
        for c in f"{sample.patient}{sample.id}{lane.id}":
            if c.isspace():
                raise ValueError(f"sample {sample.id!r}: patient, sample and lane must not contain whitespace")
    groups = await asyncio.gather(
        *(read_group(sample, lane, reference.fasta, seq_platform, seq_center) for lane in sample.lanes)
    )

    async def qc(lane: Lane) -> Dir | None:
        if skip_fastqc:
            return None
        return await fastqc(reads_1=lane.fastq_1, reads_2=lane.fastq_2, prefix=f"{sample.id}-{lane.id}", args="--quiet")

    async def align(lane: Lane, rg: str) -> File:
        return await bwa_mem(
            lane.fastq_1,
            lane.fastq_2,
            reference.bwa_index,
            prefix=f"{sample.id}-{lane.id}.sorted",
            read_group=rg,
            args="-K 100000000 -Y",
        )

    reports, bams = await asyncio.gather(
        asyncio.gather(*(qc(lane) for lane in sample.lanes)),
        asyncio.gather(*(align(lane, rg) for lane, (rg, _) in zip(sample.lanes, groups))),
    )
    return MappedSample(
        sample=sample,
        lane_bams=list(bams),
        read_groups=[rg for rg, _ in groups],
        fastqc=[r for r in reports if r is not None],
    )


@dataclass
class AlignmentQC:
    """samtools stats and mosdepth on a CRAM."""

    samtools_stats: File  # <cram name>.stats
    mosdepth: Dir  # <sample>.<stage>.mosdepth.*


@dataclass
class DuplicatesMarked:
    """The sample's alignment going into recalibration: duplicate-marked, or (skip_markduplicates) just merged."""

    cram: File  # <sample>.md.cram, or <sample>.sorted.cram when duplicate marking is skipped
    crai: File
    metrics: File | None  # <sample>.md.cram.metrics; None when duplicate marking is skipped
    qc: AlignmentQC


async def alignment_qc(
    sample_id: str, stage: str, cram: File, crai: File, reference: Reference, wes: bool
) -> AlignmentQC:
    """Upstream's CRAM QC: samtools stats, and mosdepth (500 bp windows for WGS, the intervals for WES)."""
    stats, depth = await asyncio.gather(
        samtools_stats(bam=cram, fasta=reference.fasta, fai=reference.fai),
        mosdepth(
            cram,
            crai,
            prefix=f"{sample_id}.{stage}",
            bed=reference.intervals if wes else None,
            fasta=reference.fasta,
            args="" if wes else "-n --fast-mode --by 500",
        ),
    )
    return AlignmentQC(samtools_stats=stats, mosdepth=depth)


async def mark_duplicates(mapped: MappedSample, reference: Reference, wes: bool = False) -> DuplicatesMarked:
    """MarkDuplicates across the sample's lane BAMs into ``<sample>.md.cram``, then QC it, as upstream.

    Duplicates are flagged, not removed. ``wes`` (exome / targeted data)
    restricts mosdepth to the calling intervals; for WGS it reports 500 bp
    windows.
    """
    sample_id = mapped.sample.id
    marked = await gatk4_markduplicates(
        mapped.lane_bams,
        f"{sample_id}.md.cram",
        fasta=reference.fasta,
        fai=reference.fai,
        args="--REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT",
    )
    assert marked.index is not None
    qc = await alignment_qc(sample_id, "md", marked.alignment, marked.index, reference, wes)
    return DuplicatesMarked(cram=marked.alignment, crai=marked.index, metrics=marked.metrics, qc=qc)


@dataclass
class Recalibrated:
    cram: File  # <sample>.recal.cram
    crai: File
    table: File  # the (gathered) recalibration table, <sample>.recal.table
    qc: AlignmentQC


def simple_name(f: File) -> str:
    """A file's name without any extension (Nextflow's ``simpleName``)."""
    return name(f).split(".", 1)[0]


async def recalibrate(
    marked: DuplicatesMarked, sample_id: str, reference: Reference, wes: bool = False
) -> Recalibrated:
    """Base-quality score recalibration of the duplicate-marked CRAM, as upstream.

    BaseRecalibrator runs per interval chunk against the known sites (dbSNP
    and the known indels), the tables are gathered, ApplyBQSR runs per chunk
    and the chunks are merged into ``<sample>.recal.cram``, which is then
    QC'd. Without intervals it all runs once over the whole genome.
    """
    known_sites = ([reference.dbsnp] if reference.dbsnp else []) + reference.known_indels
    known_sites_tbi = ([reference.dbsnp_tbi] if reference.dbsnp_tbi else []) + reference.known_indels_tbi
    if not known_sites:
        raise ValueError("base-quality recalibration needs known sites: give dbsnp and/or known_indels")
    chunks: list[File | None] = list(reference.intervals_split) or [None]
    single = len(chunks) == 1

    def prefix(chunk: File | None) -> str:
        # As upstream: <sample>.recal, or <sample>_<interval>.recal per chunk when scattered.
        return f"{sample_id}.recal" if single or chunk is None else f"{sample_id}_{simple_name(chunk)}.recal"

    ref = (reference.fasta, reference.fai, reference.sequence_dict)
    tables = await asyncio.gather(
        *(
            gatk4_baserecalibrator(
                marked.cram, marked.crai, *ref, known_sites, known_sites_tbi, prefix=prefix(chunk), intervals=chunk
            )
            for chunk in chunks
        )
    )
    table = tables[0] if single else await gatk4_gatherbqsrreports(list(tables), prefix=f"{sample_id}.recal")

    crams = await asyncio.gather(
        *(
            gatk4_applybqsr(marked.cram, marked.crai, table, *ref, prefix=prefix(chunk), intervals=chunk)
            for chunk in chunks
        )
    )
    cram = (
        crams[0]
        if single
        else await samtools_merge(
            list(crams),
            prefix=f"{sample_id}.recal",
            fasta=reference.fasta,
            fai=reference.fai,
            args="--output-fmt cram,version=3.0",
        )
    )
    crai = await samtools_index(bam=cram)
    qc = await alignment_qc(sample_id, "recal", cram, crai, reference, wes)
    return Recalibrated(cram=cram, crai=crai, table=table, qc=qc)


async def vcf_has_records(vcf: File) -> bool:
    """Whether a (b)gzipped VCF has any variant records, reading only as far as the first one."""
    decompressor = zlib.decompressobj(zlib.MAX_WBITS | 16)
    pending = b""
    async with vcf.open("rb") as fh:
        while True:
            chunk = bytes(await fh.read(64 * 1024))
            if not chunk:
                return False
            # bgzf is a series of gzip members: restart the decompressor at each one.
            while chunk:
                pending += decompressor.decompress(chunk)
                chunk = decompressor.unused_data
                if decompressor.eof:
                    decompressor = zlib.decompressobj(zlib.MAX_WBITS | 16)
            lines = pending.split(b"\n")
            pending = lines.pop()
            if any(line and not line.startswith(b"#") for line in lines):
                return True


@dataclass
class GermlineCalls:
    calls: IndexedVcf  # <sample>.haplotypecaller.vcf.gz (merged across interval chunks)
    scored: IndexedVcf | None  # <sample>.cnn.vcf.gz; None when the filter is skipped
    filtered: IndexedVcf | None  # <sample>.haplotypecaller.filtered.vcf.gz; None when skipped

    @property
    def final(self) -> IndexedVcf:
        """The VCF downstream steps use: the filtered one when filtering ran."""
        return self.filtered or self.calls


async def call_variants(
    cram: File,
    crai: File,
    sample_id: str,
    reference: Reference,
    skip_filter: bool = False,
    pcr_indel_model: str = "CONSERVATIVE",
    filter_args: str = "--info-key CNN_1D",
) -> GermlineCalls:
    """GATK HaplotypeCaller germline calling of a (recalibrated) CRAM, as upstream.

    HaplotypeCaller runs per interval chunk (annotating dbSNP IDs), the
    chunks are merged into ``<sample>.haplotypecaller.vcf.gz``, and unless
    ``skip_filter`` the calls are scored with GATK's 1D CNN over the calling
    intervals and filtered by sensitivity tranches against dbSNP and the
    known indels (``filter_args``; upstream's test profile adds
    ``--indel-tranche 0``). GATK refuses to score or filter a VCF without
    records, which upstream lets fail; a sample with no calls skips the
    filter instead (``scored`` and ``filtered`` are None).
    """
    chunks: list[File | None] = list(reference.intervals_split) or [None]
    single = len(chunks) == 1
    ref = (reference.fasta, reference.fai, reference.sequence_dict)

    def prefix(chunk: File | None) -> str:
        # As upstream: <sample>.haplotypecaller, or .haplotypecaller.<interval> per chunk when scattered.
        if single or chunk is None:
            return f"{sample_id}.haplotypecaller"
        return f"{sample_id}.haplotypecaller.{base_name(chunk)}"

    per_chunk = await asyncio.gather(
        *(
            gatk4_haplotypecaller(
                cram,
                crai,
                *ref,
                prefix=prefix(chunk),
                intervals=chunk,
                dbsnp=reference.dbsnp,
                dbsnp_tbi=reference.dbsnp_tbi,
                args=f"--pcr-indel-model {pcr_indel_model}" if pcr_indel_model else "",
            )
            for chunk in chunks
        )
    )
    calls = (
        per_chunk[0]
        if single
        else await gatk4_mergevcfs(
            [c.vcf for c in per_chunk], prefix=f"{sample_id}.haplotypecaller", dict=reference.sequence_dict
        )
    )
    if skip_filter or not await vcf_has_records(calls.vcf):
        return GermlineCalls(calls=calls, scored=None, filtered=None)

    known_sites = ([reference.dbsnp] if reference.dbsnp else []) + reference.known_indels
    known_sites_tbi = ([reference.dbsnp_tbi] if reference.dbsnp_tbi else []) + reference.known_indels_tbi
    if not known_sites:
        raise ValueError("HaplotypeCaller filtering needs known sites: give dbsnp and/or known_indels, or skip_filter")
    # As upstream: one pass over all calling intervals (scattering fails on chunks without SNPs).
    scored = await gatk4_cnnscorevariants(calls.vcf, calls.tbi, *ref, prefix=sample_id, intervals=reference.intervals)
    filtered = await gatk4_filtervarianttranches(
        scored.vcf,
        scored.tbi,
        known_sites,
        known_sites_tbi,
        *ref,
        prefix=f"{sample_id}.haplotypecaller",
        args=filter_args,
    )
    return GermlineCalls(calls=calls, scored=scored, filtered=filtered)


@dataclass
class VcfQC:
    bcftools_stats: File  # <vcf name>.bcftools_stats.txt
    tstv_count: File  # <vcf name>.TsTv.count
    tstv_qual: File  # <vcf name>.TsTv.qual
    filter_summary: File  # <vcf name>.FILTER.summary


def vcf_prefix(vcf: File) -> str:
    """A VCF's name without ``.vcf[.gz]`` (upstream's ``vcf.baseName - ".vcf"``)."""
    n = name(vcf)
    return n.removesuffix(".gz").removesuffix(".vcf")


async def vcf_qc(vcf: IndexedVcf) -> VcfQC:
    """QC of a sample's calls, as upstream: bcftools stats and vcftools Ts/Tv and FILTER summaries."""
    prefix = vcf_prefix(vcf.vcf)
    stats, count, qual, summary = await asyncio.gather(
        bcftools_stats(vcf.vcf, prefix=prefix),
        vcftools(vcf.vcf, prefix=prefix, args="--TsTv-by-count"),
        vcftools(vcf.vcf, prefix=prefix, args="--TsTv-by-qual"),
        vcftools(vcf.vcf, prefix=prefix, args="--FILTER-summary"),
    )

    async def pick(results: Dir, suffix: str) -> File:
        f = await results.get_file(f"{prefix}.{suffix}")
        if f is None:
            raise FileNotFoundError(f"vcftools wrote no {prefix}.{suffix}")
        return f

    return VcfQC(
        bcftools_stats=stats,
        tstv_count=await pick(count, "TsTv.count"),
        tstv_qual=await pick(qual, "TsTv.qual"),
        filter_summary=await pick(summary, "FILTER.summary"),
    )


@dataclass
class SampleQC:
    """What one sample contributes to the MultiQC report."""

    sample_id: str
    fastqc: list[Dir] = field(default_factory=list)
    markduplicates: DuplicatesMarked | None = None
    recalibrated: Recalibrated | None = None
    vcf: VcfQC | None = None


async def multiqc_report(samples: list[SampleQC]) -> MultiQCResult:
    """One MultiQC report over every sample's QC, with upstream's config.

    The reports are laid out under the names upstream gives them, which is
    what MultiQC's modules and the config's sample-name cleaning key on.
    """
    root = Path(tempfile.mkdtemp(prefix="multiqc_"))
    downloads = []

    def put(f: File | None, rel: str) -> None:
        if f is not None:
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            downloads.append(f.download(str(root / rel)))

    async def put_dir(d: Dir | None, rel: str, suffixes: tuple[str, ...] = ()) -> None:
        if d is None:
            return
        base = d.path.rstrip("/")
        async for f in d.walk():
            sub = f.path[len(base) :].lstrip("/")
            if not suffixes or sub.endswith(suffixes):
                put(f, f"{rel}/{sub}")

    for s in samples:
        sid = s.sample_id
        for i, report in enumerate(s.fastqc):
            await put_dir(report, f"fastqc/{sid}/{i}", ("_fastqc.zip",))
        md = s.markduplicates
        if md is not None:
            put(md.metrics, f"markduplicates/{name(md.metrics) if md.metrics else ''}")
            put(md.qc.samtools_stats, f"samtools/{name(md.cram)}.stats")
            await put_dir(md.qc.mosdepth, f"mosdepth/{sid}", (".txt",))
        recal = s.recalibrated
        if recal is not None:
            put(recal.table, f"bqsr/{name(recal.table)}")
            put(recal.qc.samtools_stats, f"samtools/{sid}.recal.cram.stats")
            await put_dir(recal.qc.mosdepth, f"mosdepth/{sid}", (".txt",))
        q = s.vcf
        if q is not None:
            for f in (q.bcftools_stats, q.tstv_count, q.tstv_qual, q.filter_summary):
                put(f, f"vcf/{sid}/{name(f)}")

    await asyncio.gather(*downloads)
    data = await Dir.from_local(str(root))
    config = await File.from_local(str(asset("multiqc_config_variant_calling.yml")))
    return await multiqc(data, configs=[config], prefix="multiqc_report")


async def merge_lanes(mapped: MappedSample, reference: Reference, wes: bool = False) -> DuplicatesMarked:
    """Without duplicate marking, as upstream: the lane BAMs merged into ``<sample>.sorted.cram``, then QC'd."""
    sample_id = mapped.sample.id
    cram = await samtools_merge(
        mapped.lane_bams,
        prefix=f"{sample_id}.sorted",
        fasta=reference.fasta,
        fai=reference.fai,
        args="--output-fmt cram,version=3.0",
    )
    crai = await samtools_index(bam=cram)
    qc = await alignment_qc(sample_id, "sorted", cram, crai, reference, wes)
    return DuplicatesMarked(cram=cram, crai=crai, metrics=None, qc=qc)


@dataclass
class SampleResult:
    sample: Sample
    mapped: MappedSample
    markduplicates: DuplicatesMarked  # merged lanes only, when skip_markduplicates
    recalibrated: Recalibrated | None  # None when skip_baserecalibrator
    calls: GermlineCalls
    vcf_qc: VcfQC | None  # None when skip_vcf_qc


@dataclass
class VariantCallingResult:
    reference: Reference
    samples: list[SampleResult]
    multiqc: MultiQCResult | None = None  # None when skip_multiqc
    outdir: Dir | None = None  # the published results tree, when outdir/publish_results was asked for


async def variant_calling(
    samplesheet: File,
    fasta: File,
    fai: File | None = None,
    sequence_dict: File | None = None,
    bwa_index: Dir | None = None,
    dbsnp: File | None = None,
    dbsnp_tbi: File | None = None,
    known_indels: list[File] | None = None,
    known_indels_tbi: list[File] | None = None,
    intervals: File | None = None,
    no_intervals: bool = False,
    nucleotides_per_second: int = 200000,
    wes: bool = False,
    seq_platform: str = "ILLUMINA",
    seq_center: str = "",
    skip_fastqc: bool = False,
    skip_markduplicates: bool = False,
    skip_baserecalibrator: bool = False,
    skip_haplotypecaller_filter: bool = False,
    haplotypecaller_filter_args: str = "--info-key CNN_1D",
    gatk_pcr_indel_model: str = "CONSERVATIVE",
    skip_vcf_qc: bool = False,
    skip_multiqc: bool = False,
    outdir: str | None = None,
    publish_results: bool = False,
) -> VariantCallingResult:
    """Germline short-variant calling with GATK HaplotypeCaller, from FASTQ, as upstream's germline route.

    ``samplesheet`` is upstream's ``--input``: a CSV with columns ``patient``,
    ``sex``, ``status`` (0 = normal; tumor samples aren't supported yet),
    ``sample``, ``lane``, ``fastq_1`` and ``fastq_2`` (empty for single-end),
    one row per lane. It is read before anything runs. The reference is
    prepared once (see :func:`prepare_reference`; anything not supplied is
    built), and every sample is mapped, duplicate-marked, recalibrated
    against the known sites, called and QC'd in parallel, then reported
    together in MultiQC.

    ``intervals`` (BED or interval list) restricts calling to target regions
    and sets the scatter-gather chunks (sized by ``nucleotides_per_second``);
    ``wes`` marks exome / targeted data (mosdepth reports on the intervals
    rather than 500 bp windows). Base-quality recalibration and the CNN
    filter need known sites (``dbsnp`` and/or ``known_indels``) unless
    skipped.

    Outputs: everything is returned in the :class:`VariantCallingResult`.
    Like upstream's ``--outdir``, ``outdir`` (e.g. ``s3://my-bucket/vc``)
    also lays the results out there in upstream's folder structure
    (:func:`variant_calling_layout`); ``publish_results=True`` does the same
    in Flyte's own storage. Either way ``result.outdir`` is that tree as a
    Dir, for downstream tasks.
    """
    samples = samples_from_samplesheet(await read_samplesheet(samplesheet))
    if not samples:
        raise ValueError("the samplesheet has no samples")
    if not skip_baserecalibrator and dbsnp is None and not known_indels:
        raise ValueError("base-quality recalibration needs dbsnp and/or known_indels (or skip_baserecalibrator)")
    if not skip_haplotypecaller_filter and dbsnp is None and not known_indels:
        raise ValueError("the HaplotypeCaller filter needs dbsnp and/or known_indels (or skip_haplotypecaller_filter)")

    reference = await prepare_reference(
        fasta,
        fai=fai,
        sequence_dict=sequence_dict,
        bwa_index=bwa_index,
        dbsnp=dbsnp,
        dbsnp_tbi=dbsnp_tbi,
        known_indels=known_indels,
        known_indels_tbi=known_indels_tbi,
        intervals=intervals,
        no_intervals=no_intervals,
        nucleotides_per_second=nucleotides_per_second,
    )

    async def run_sample(sample: Sample) -> SampleResult:
        mapped = await map_reads(sample, reference, seq_platform, seq_center, skip_fastqc=skip_fastqc)
        if skip_markduplicates:
            marked = await merge_lanes(mapped, reference, wes)
        else:
            marked = await mark_duplicates(mapped, reference, wes)
        recal = None if skip_baserecalibrator else await recalibrate(marked, sample.id, reference, wes)
        cram, crai = (recal.cram, recal.crai) if recal else (marked.cram, marked.crai)
        calls = await call_variants(
            cram,
            crai,
            sample.id,
            reference,
            skip_filter=skip_haplotypecaller_filter,
            pcr_indel_model=gatk_pcr_indel_model,
            filter_args=haplotypecaller_filter_args,
        )
        qc = None if skip_vcf_qc else await vcf_qc(calls.final)
        return SampleResult(sample, mapped, marked, recal, calls, qc)

    results = await asyncio.gather(*(run_sample(s) for s in samples))
    result = VariantCallingResult(reference=reference, samples=list(results))
    if not skip_multiqc:
        result.multiqc = await multiqc_report(
            [SampleQC(r.sample.id, r.mapped.fastqc, r.markduplicates, r.recalibrated, r.vcf_qc) for r in results]
        )
    if outdir is not None or publish_results:
        result.outdir = await publish(await variant_calling_layout(result), outdir)
    return result


async def variant_calling_layout(result: VariantCallingResult) -> dict[str, File | Dir | None]:
    """Upstream's results layout for a run: relative path -> output, for :func:`~flyte_bio.publish.publish`.

    Follows upstream's default ``publishDir`` paths: ``preprocessing/``
    (``markduplicates/``, ``recal_table/``, ``recalibrated/`` per sample),
    ``reports/`` (``fastqc/``, ``markduplicates/``, ``samtools/``,
    ``mosdepth/``, ``bcftools/`` and ``vcftools/`` per caller and sample),
    ``variant_calling/haplotypecaller/<sample>/`` and ``multiqc/``. Like
    upstream by default, it leaves out the lane BAMs and the prepared
    reference.
    """
    out: dict[str, File | Dir | None] = {}
    for r in result.samples:
        sid = r.sample.id
        for lane, report in zip(r.sample.lanes, r.mapped.fastqc):
            out.update(await expand(f"reports/fastqc/{sid}-{lane.id}", report))
        md = r.markduplicates
        stage = "markduplicates" if md.metrics is not None else "mapped"
        out[f"preprocessing/{stage}/{sid}/{name(md.cram)}"] = md.cram
        out[f"preprocessing/{stage}/{sid}/{name(md.cram)}.crai"] = md.crai
        if md.metrics is not None:
            out[f"reports/markduplicates/{sid}/{name(md.metrics)}"] = md.metrics
        out[f"reports/samtools/{sid}/{name(md.cram)}.stats"] = md.qc.samtools_stats
        out.update(await expand(f"reports/mosdepth/{sid}", md.qc.mosdepth))
        if r.recalibrated is not None:
            rc = r.recalibrated
            out[f"preprocessing/recal_table/{sid}/{name(rc.table)}"] = rc.table
            out[f"preprocessing/recalibrated/{sid}/{name(rc.cram)}"] = rc.cram
            out[f"preprocessing/recalibrated/{sid}/{name(rc.cram)}.crai"] = rc.crai
            out[f"reports/samtools/{sid}/{name(rc.cram)}.stats"] = rc.qc.samtools_stats
            out.update(await expand(f"reports/mosdepth/{sid}", rc.qc.mosdepth))
        calls_dir = f"variant_calling/haplotypecaller/{sid}"
        for vcf in (r.calls.calls, r.calls.filtered):
            if vcf is not None:
                out[f"{calls_dir}/{name(vcf.vcf)}"] = vcf.vcf
                out[f"{calls_dir}/{name(vcf.tbi)}"] = vcf.tbi
        if r.vcf_qc is not None:
            q = r.vcf_qc
            out[f"reports/bcftools/haplotypecaller/{sid}/{name(q.bcftools_stats)}"] = q.bcftools_stats
            for f in (q.tstv_count, q.tstv_qual, q.filter_summary):
                out[f"reports/vcftools/haplotypecaller/{sid}/{name(f)}"] = f
    if result.multiqc is not None:
        out.update(await expand("multiqc", result.multiqc.results))
    return out
