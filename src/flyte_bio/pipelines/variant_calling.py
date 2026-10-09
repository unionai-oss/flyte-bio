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

Upstream splits FASTQs into 50M-read chunks with fastp by default to map
them in parallel; that isn't ported (each lane maps as one job), which
changes how the work is spread, not the alignments.

These are plain async functions, not tasks: they run inside the caller's
task and fan out to the module tasks, so the caller's
:class:`flyte.TaskEnvironment` must ``depends_on`` :data:`flyte_bio.modules.env`.
"""

import asyncio
import zlib
from dataclasses import dataclass, field

from flyte.io import Dir, File

from flyte_bio.modules.bwa import bwa_index as build_bwa_index
from flyte_bio.modules.bwa import bwa_mem
from flyte_bio.modules.fastqc import fastqc
from flyte_bio.modules.gatk4 import gatk4_createsequencedictionary, gatk4_intervallisttobed, gatk4_markduplicates
from flyte_bio.modules.htslib import htslib_bgziptabix
from flyte_bio.modules.intervals import build_intervals, create_intervals_bed
from flyte_bio.modules.mosdepth import mosdepth
from flyte_bio.modules.samtools import samtools_faidx, samtools_stats
from flyte_bio.publish import name
from flyte_bio.samplesheet import Row


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
    """The first line of a gzipped FASTQ, reading only as much as it takes."""
    decompressor = zlib.decompressobj(zlib.MAX_WBITS | 16)
    text = b""
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
    cram: File  # <sample>.md.cram
    crai: File
    metrics: File  # <sample>.md.cram.metrics
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
