"""rnaseq — RNA-seq quantification (STAR alignment + salmon), rnaseq 3.26.0.

Composes :mod:`flyte_bio.modules` tasks into the ``star_salmon`` path:

- :func:`prepare_genome` — decompress references, filter the GTF to the
  genome's sequences, append an additional FASTA, derive the gene BED,
  transcript FASTA and chrom sizes, and build the STAR index.
- :func:`preprocess_reads` — merge a sample's runs, lint them, FastQC the
  raw reads, trim with Trim Galore (FastQC on the trimmed reads), lint again
  and count the reads that survive (samples below ``min_trimmed_reads`` are
  dropped), then optionally BBSplit away reads from other genomes (and lint
  the result) — all as upstream.
- :func:`infer_strandedness` — for ``strandedness="auto"`` samples, subsample
  the trimmed reads and run salmon (``--skipQuant``) to classify the library
  as forward / reverse / unstranded, as upstream does.
- :func:`align_star` — align the trimmed reads with STAR (emitting a
  transcriptome BAM), then sort / index / stats the genome BAM.
- :func:`mark_duplicates` — Picard MarkDuplicates on the genome BAM, then
  index / stats the marked BAM, which replaces it downstream (as upstream).
- :func:`bigwig_coverage` — genome coverage bigWigs (per strand for stranded
  libraries, plus combined), as upstream.
- StringTie (reference-guided, ``-e``), dupRadar, Qualimap (on a name-sorted
  copy), RSeQC and the featureCounts biotype QC on the (marked) genome BAM.
- :func:`quantify_salmon_bam` — salmon alignment-mode quantification of the
  transcriptome BAM.
- :func:`merge_quantifications` — tx2gene + tximport across all samples into
  gene/transcript count, TPM and length matrices, bundled as gene- and
  transcript-level SummarizedExperiment RDS files, then DESeq2 PCA /
  sample-distance QC on the merged gene counts.
- :func:`rnaseq` — run all of the above, samples in parallel.

These are plain async functions, not tasks: they run inside the caller's
task and fan out to the module tasks, so the caller's
:class:`flyte.TaskEnvironment` must ``depends_on`` :data:`flyte_bio.modules.env`.

- :func:`multiqc_report` — one MultiQC report over every QC output, with
  upstream's config and custom content (failed-sample tables, strandedness
  checks, paired-end sample merging, FASTQ-to-sample name replacements).
  Nextflow-specific sections (run parameters, software versions, methods
  text) are omitted. rRNA removal and UMI deduplication (off upstream by
default) aren't ported either.
"""

import asyncio
import csv
import json
import math
import re
import tempfile
from dataclasses import dataclass, field, replace
from pathlib import Path

import yaml
from flyte.io import Dir, File

from flyte_bio.modules.bbmap import BBSplitResult, bbsplit, bbsplit_index
from flyte_bio.modules.bedtools import bedtools_genomecov
from flyte_bio.modules.cat import cat_fastq
from flyte_bio.modules.catadditionalfasta import cat_additional_fasta
from flyte_bio.modules.deseq2_qc import Deseq2QCResult, deseq2_qc
from flyte_bio.modules.dupradar import DupradarResult, dupradar
from flyte_bio.modules.fastqc import fastqc
from flyte_bio.modules.fq import fq_lint, fq_subsample
from flyte_bio.modules.gffread import gffread_gff_to_gtf, gffread_transcripts_fasta
from flyte_bio.modules.gtf2bed import gtf2bed
from flyte_bio.modules.gtffilter import gtf_filter
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.multiqc import MultiQCResult, multiqc
from flyte_bio.modules.multiqccustombiotype import BiotypeQC, multiqc_custom_biotype
from flyte_bio.modules.picard import picard_markduplicates
from flyte_bio.modules.qualimap import qualimap_rnaseq
from flyte_bio.modules.rseqc import DEFAULT_MODULES as RSEQC_MODULES
from flyte_bio.modules.rseqc import rseqc
from flyte_bio.modules.salmon import salmon_index, salmon_quant_bam, salmon_quant_reads
from flyte_bio.modules.samtools import (
    samtools_faidx,
    samtools_flagstat,
    samtools_idxstats,
    samtools_index,
    samtools_sort,
    samtools_stats,
)
from flyte_bio.modules.star import star_align, star_genome_generate
from flyte_bio.modules.stringtie import StringTieResult, stringtie
from flyte_bio.modules.subread import FeatureCountsResult, featurecounts
from flyte_bio.modules.summarizedexperiment import summarized_experiment
from flyte_bio.modules.trimgalore import TrimGaloreResult, trimgalore
from flyte_bio.modules.tx2gene import tx2gene
from flyte_bio.modules.tximport import TximportResult, collect_quants, tximport
from flyte_bio.modules.ucsc import bedclip, bedgraphtobigwig
from flyte_bio.scripts import path as asset

STRANDEDNESS = ("auto", "forward", "reverse", "unstranded")

# Upstream's BBSplit options (index build and split alike).
BBSPLIT_ARGS = "build=1 ambiguous2=all maxindel=150000 ow=f"


@dataclass
class Sample:
    """One sample: its sequencing runs (merged before alignment) and strandedness.

    ``fastq_2`` is empty for single-end samples; otherwise it pairs
    one-to-one with ``fastq_1`` in the same run order.
    """

    id: str
    fastq_1: list[File]
    fastq_2: list[File] = field(default_factory=list)
    strandedness: str = "auto"

    def __post_init__(self) -> None:
        if not self.fastq_1:
            raise ValueError(f"sample {self.id!r} has no reads")
        if self.fastq_2 and len(self.fastq_2) != len(self.fastq_1):
            raise ValueError(
                f"sample {self.id!r} mixes single- and paired-end runs "
                f"({len(self.fastq_1)} R1 vs {len(self.fastq_2)} R2)"
            )
        if self.strandedness not in STRANDEDNESS:
            raise ValueError(f"sample {self.id!r}: strandedness must be one of {STRANDEDNESS}")

    @property
    def single_end(self) -> bool:
        return not self.fastq_2


@dataclass
class Genome:
    fasta: File
    gtf: File
    fai: File
    chrom_sizes: File
    gene_bed: File
    transcript_fasta: File
    star_index: Dir


@dataclass
class StarAlignment:
    bam: File  # coordinate-sorted genome BAM
    bai: File
    stats: File
    flagstat: File
    idxstats: File
    transcriptome_bam: File
    log_final: File
    splice_junctions: File
    star_output: Dir  # everything STAR wrote, for logs not surfaced above


@dataclass
class PreprocessedReads:
    reads_1: File  # reads to align: after trimming / BBSplit, whichever ran last
    reads_2: File | None
    raw_fastqc: Dir | None  # <sample>_raw[_1|_2]_fastqc.{html,zip}
    trimming: TrimGaloreResult | None  # reports + trimmed-read FastQC
    reads_after_trimming: float | None
    bbsplit: BBSplitResult | None = None
    lint: dict[str, File] = field(default_factory=dict)  # stage (raw/trimmed/bbsplit) -> fq lint log


@dataclass
class StrandednessAnalysis:
    """Upstream's salmon-based strandedness call for an ``auto`` sample."""

    inferred: str  # forward / reverse / unstranded / undetermined
    forward_pct: float  # % of fragments
    reverse_pct: float
    unstranded_pct: float
    lib_format_counts: File  # salmon's lib_format_counts.json for the subsample


@dataclass
class MarkedDuplicates:
    bam: File  # <sample>.markdup.sorted.bam (duplicates flagged, not removed)
    bai: File
    metrics: File  # Picard duplication metrics
    stats: File
    flagstat: File
    idxstats: File


@dataclass
class SampleResult:
    sample: str
    percent_mapped: float  # STAR "Uniquely mapped reads %"
    mapping_passed: bool  # percent_mapped >= min_mapped_reads; if not, no genome-BAM steps ran
    strandedness: str  # as used for quantification (inferred for ``auto`` samples)
    strandedness_analysis: StrandednessAnalysis | None  # set for ``auto`` samples
    preprocessing: PreprocessedReads
    alignment: StarAlignment
    salmon: Dir
    markduplicates: MarkedDuplicates | None = None  # None when skip_markduplicates
    dupradar: DupradarResult | None = None  # None when skip_dupradar
    qualimap: Dir | None = None  # Qualimap rnaseq report dir; None when skip_qualimap
    rseqc: dict[str, Dir] = field(default_factory=dict)  # RSeQC module -> outputs; empty when skip_rseqc
    biotype_counts: FeatureCountsResult | None = None  # featureCounts by biotype; None when skip_biotype_qc
    biotype_qc: BiotypeQC | None = None  # MultiQC biotype tables
    stringtie: StringTieResult | None = None  # None when skip_stringtie
    bigwig: dict[str, File] = field(default_factory=dict)  # forward/reverse/combined -> bigWig; empty when skip_bigwig


@dataclass
class MergedQuantification:
    tx2gene: File
    tximport: TximportResult
    gene_rds: File  # gene-level SummarizedExperiment
    transcript_rds: File  # transcript-level SummarizedExperiment


@dataclass
class RnaseqResult:
    genome: Genome
    samples: list[SampleResult]
    salmon: MergedQuantification
    # Samples dropped for too few reads after trimming -> surviving read count.
    failed_trimming: dict[str, float] = field(default_factory=dict)
    deseq2_qc: Deseq2QCResult | None = None  # None when skip_deseq2_qc
    multiqc: MultiQCResult | None = None  # None when skip_multiqc


async def maybe_gunzip(file: File | None) -> File | None:
    """Decompress ``file`` when its name ends in ``.gz`` (as upstream decides)."""
    if file is not None and file.path.endswith(".gz"):
        return await gunzip(archive=file)
    return file


async def prepare_genome(
    fasta: File,
    gtf: File | None = None,
    gff: File | None = None,
    transcript_fasta: File | None = None,
    additional_fasta: File | None = None,
    biotype: str = "gene_biotype",
) -> Genome:
    """Build the reference set the ``star_salmon`` path needs.

    Supply ``gtf`` or ``gff`` (a GFF is converted). Without a
    ``transcript_fasta`` one is extracted from the genome with gffread
    (upstream's default here is RSEM's extractor, not yet ported).
    ``biotype`` is the GTF attribute written for ``additional_fasta``
    features.
    """
    if gtf is None and gff is None:
        raise ValueError("prepare_genome needs a gtf or a gff")

    fasta, gtf, gff, transcript_fasta, additional_fasta = await asyncio.gather(
        maybe_gunzip(fasta),
        maybe_gunzip(gtf),
        maybe_gunzip(gff),
        maybe_gunzip(transcript_fasta),
        maybe_gunzip(additional_fasta),
    )
    if gtf is None:
        gtf = await gffread_gff_to_gtf(gff=gff)

    gtf = await gtf_filter(gtf=gtf, fasta=fasta)
    if additional_fasta is not None:
        fasta, gtf = await cat_additional_fasta(fasta=fasta, gtf=gtf, add_fasta=additional_fasta, biotype=biotype)

    async def transcripts() -> File:
        if transcript_fasta is not None:
            return transcript_fasta
        return await gffread_transcripts_fasta(annotation=gtf, fasta=fasta)

    gene_bed, transcripts_out, (fai, chrom_sizes), star_index = await asyncio.gather(
        gtf2bed(gtf=gtf),
        transcripts(),
        samtools_faidx(fasta=fasta),
        star_genome_generate(fasta=fasta, gtf=gtf),
    )
    return Genome(
        fasta=fasta,
        gtf=gtf,
        fai=fai,
        chrom_sizes=chrom_sizes,
        gene_bed=gene_bed,
        transcript_fasta=transcripts_out,
        star_index=star_index,
    )


def star_align_args(sample: Sample, seq_platform: str = "", seq_center: str = "") -> str:
    """Upstream's STAR arguments for the ``star_salmon`` aligner.

    The read-group fields are space-separated and unquoted: the shell task
    passes args through word splitting, so quotes would end up literally in
    the BAM header. Sample IDs therefore must not contain whitespace.
    """
    rg = [f"ID:{sample.id}", f"SM:{sample.id}"]
    if seq_platform:
        rg.append(f"PL:{seq_platform}")
    if seq_center:
        rg.append(f"CN:{seq_center}")
    return " ".join(
        [
            "--quantMode TranscriptomeSAM",
            "--outSAMtype BAM Unsorted",
            "--outSAMattributes NH HI AS NM MD",
            "--readFilesCommand zcat",
            "--twopassMode Basic",
            "--runRNGseed 0",
            "--outFilterMultimapNmax 20",
            "--alignSJDBoverhangMin 1",
            "--outSAMstrandField intronMotif",
            "--quantTranscriptomeSAMoutput BanSingleEnd",
            "--outSAMattrRGline " + " ".join(rg),
        ]
    )


async def star_file(output: Dir, name: str) -> File:
    found = await output.get_file(name)
    if found is None:
        raise FileNotFoundError(f"STAR output {output.path} has no {name}")
    return found


async def lint(reads_1: File, reads_2: File | None, args: str = "") -> File:
    """fq lint the reads; the task fails (failing the run) on invalid input."""
    return await fq_lint(reads_1=reads_1, reads_2=reads_2, args=args)


async def preprocess_reads(
    sample: Sample,
    opts: "RunOptions",
    bbsplit_index_dir: "asyncio.Task[Dir] | None" = None,
) -> PreprocessedReads:
    """Merge, lint, QC, trim, filter and (optionally) BBSplit one sample's reads, as upstream.

    When trimming leaves fewer than ``opts.min_trimmed_reads`` reads the
    sample stops there (BBSplit is skipped); the caller drops it.
    """
    skip_fastqc, skip_trimming = opts.skip_fastqc, opts.skip_trimming
    # Upstream's merged-read names: <sample>.merged.fastq.gz / <sample>_{1,2}.merged.fastq.gz
    mate_1, mate_2 = ("", "") if sample.single_end else ("_1", "_2")
    reads_1, reads_2 = await asyncio.gather(
        cat_fastq(sample.fastq_1, prefix=f"{sample.id}{mate_1}.merged"),
        cat_fastq(sample.fastq_2, prefix=f"{sample.id}{mate_2}.merged")
        if sample.fastq_2
        else asyncio.sleep(0, result=None),
    )
    # Upstream's `withName: 'FQ_LINT'` args also reach its aliased lint steps
    # (after trimming / BBSplit), so every lint gets extra_fqlint_args — by
    # default `--disable-validator P001`, since SRA-style headers carry /1 /2.
    lint_logs: dict[str, File] = {}
    if not opts.skip_linting:
        lint_logs["raw"] = await lint(reads_1, reads_2, opts.extra_fqlint_args)

    async def raw_qc() -> Dir | None:
        if skip_fastqc:
            return None
        return await fastqc(
            reads_1=reads_1,
            reads_2=reads_2,
            prefix=f"{sample.id}_raw",
            args="--quiet",
        )

    async def trim() -> TrimGaloreResult | None:
        if skip_trimming:
            return None
        return await trimgalore(reads_1, reads_2, prefix=f"{sample.id}_trimmed", fastqc=not skip_fastqc)

    raw_fastqc, trimming = await asyncio.gather(raw_qc(), trim())
    reads_after_trimming = None
    if trimming is not None:
        reads_1, reads_2 = trimming.reads_1, trimming.reads_2
        if not opts.skip_linting:
            lint_logs["trimmed"] = await lint(reads_1, reads_2, opts.extra_fqlint_args)
        reads_after_trimming = await trimming.reads_after_filtering()
    pre = PreprocessedReads(reads_1, reads_2, raw_fastqc, trimming, reads_after_trimming, lint=lint_logs)
    if reads_after_trimming is not None and reads_after_trimming < opts.min_trimmed_reads:
        return pre

    if bbsplit_index_dir is not None:
        split = await bbsplit(reads_1, reads_2, await bbsplit_index_dir, prefix=sample.id, args=BBSPLIT_ARGS)
        pre.bbsplit, pre.reads_1, pre.reads_2 = split, split.reads_1, split.reads_2
        if not opts.skip_linting:
            lint_logs["bbsplit"] = await lint(split.reads_1, split.reads_2, opts.extra_fqlint_args)
    return pre


async def align_star(
    sample: Sample,
    reads_1: File,
    reads_2: File | None,
    genome: Genome,
    seq_platform: str = "",
    seq_center: str = "",
) -> StarAlignment:
    """Align one sample's (trimmed) reads with STAR, then sort / index / stats the BAM."""
    if any(c.isspace() for c in sample.id):
        raise ValueError(f"sample id {sample.id!r} must not contain whitespace")

    output = await star_align(
        reads_1=reads_1,
        reads_2=reads_2,
        index=genome.star_index,
        gtf=genome.gtf,
        args=star_align_args(sample, seq_platform, seq_center),
    )
    bam, transcriptome_bam, log_final, splice_junctions = await asyncio.gather(
        star_file(output, "Aligned.out.bam"),
        star_file(output, "Aligned.toTranscriptome.out.bam"),
        star_file(output, "Log.final.out"),
        star_file(output, "SJ.out.tab"),
    )

    bam_sorted = await samtools_sort(bam=bam)
    bai = await samtools_index(bam=bam_sorted)
    stats, flagstat, idxstats = await asyncio.gather(
        samtools_stats(bam=bam_sorted),
        samtools_flagstat(bam=bam_sorted),
        samtools_idxstats(bam=bam_sorted, bai=bai),
    )
    return StarAlignment(
        bam=bam_sorted,
        bai=bai,
        stats=stats,
        flagstat=flagstat,
        idxstats=idxstats,
        transcriptome_bam=transcriptome_bam,
        log_final=log_final,
        splice_junctions=splice_junctions,
        star_output=output,
    )


async def star_percent_mapped(log_final: File) -> float:
    """STAR's uniquely-mapped percentage from Log.final.out, parsed as upstream does (0 if absent)."""
    async with log_final.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    percent = 0.0
    for line in text.splitlines():
        if m := re.search(r"Uniquely mapped reads %\s*\|\s*([\d\.]+)%", line):
            percent = float(m.group(1))
    return percent


async def bigwig_coverage(sample: Sample, bam: File, chrom_sizes: File) -> dict[str, File]:
    """Genome coverage bigWigs, as upstream: combined, plus per strand when stranded.

    Upstream names the per-strand tracks by transcript strand, so for a
    reverse-stranded library ``<sample>.forward.bigWig`` comes from the
    ``-strand -`` reads (and the intermediate bedGraph prefixes swap too).
    """

    async def track(cov_prefix: str, cov_args: str, clip_prefix: str, bw_prefix: str) -> File:
        bedgraph = await bedtools_genomecov(bam, prefix=cov_prefix, extension="bedGraph", sort=True, args=cov_args)
        clipped = await bedclip(bedgraph, chrom_sizes, prefix=clip_prefix)
        return await bedgraphtobigwig(clipped, chrom_sizes, prefix=bw_prefix)

    rev = sample.strandedness == "reverse"
    tracks = {"combined": track(sample.id, "-split -bg", f"{sample.id}.clip", sample.id)}
    if sample.strandedness in ("forward", "reverse"):
        tracks["forward"] = track(
            f"{sample.id}.reverse" if rev else f"{sample.id}.forward",
            f"-split -du -strand {'-' if rev else '+'} -bg",
            f"{sample.id}.clip.forward",
            f"{sample.id}.forward",
        )
        tracks["reverse"] = track(
            f"{sample.id}.forward" if rev else f"{sample.id}.reverse",
            f"-split -du -strand {'+' if rev else '-'} -bg",
            f"{sample.id}.clip.reverse",
            f"{sample.id}.reverse",
        )
    results = await asyncio.gather(*tracks.values())
    return dict(zip(tracks, results))


# Upstream's Picard MarkDuplicates options for the genome BAM.
MARKDUPLICATES_ARGS = "--ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp"


async def mark_duplicates(sample: Sample, alignment: StarAlignment, genome: Genome) -> MarkedDuplicates:
    """Picard MarkDuplicates on the sorted genome BAM, then index / stats it, as upstream."""
    marked = await picard_markduplicates(
        alignment.bam,
        prefix=f"{sample.id}.markdup.sorted",
        fasta=genome.fasta,
        fai=genome.fai,
        args=MARKDUPLICATES_ARGS,
    )
    bai = await samtools_index(bam=marked.bam)
    stats, flagstat, idxstats = await asyncio.gather(
        samtools_stats(bam=marked.bam),
        samtools_flagstat(bam=marked.bam),
        samtools_idxstats(bam=marked.bam, bai=bai),
    )
    return MarkedDuplicates(marked.bam, bai, marked.metrics, stats, flagstat, idxstats)


def calculate_strandedness(
    forward: float,
    reverse: float,
    unstranded: float,
    stranded_threshold: float = 0.8,
    unstranded_threshold: float = 0.1,
) -> tuple[str, float, float, float]:
    """Upstream's strandedness call from fragment counts.

    Returns ``(strandedness, forward %, reverse %, unstranded %)``. A library
    is forward/reverse when that strand holds at least ``stranded_threshold``
    of the stranded fragments, unstranded when the two strands differ by at
    most ``unstranded_threshold``, and ``undetermined`` otherwise.
    """
    total = forward + reverse + unstranded
    stranded = forward + reverse
    strandedness = "undetermined"
    if stranded > 0:
        forward_prop, reverse_prop = forward / stranded, reverse / stranded
        if forward_prop >= stranded_threshold:
            strandedness = "forward"
        elif reverse_prop >= stranded_threshold:
            strandedness = "reverse"
        elif abs(forward_prop - reverse_prop) <= unstranded_threshold:
            strandedness = "unstranded"
    pct = (lambda n: n / total * 100) if total else (lambda n: 0.0)
    return strandedness, pct(forward), pct(reverse), pct(unstranded)


def salmon_strandedness(
    lib_format_counts: dict, stranded_threshold: float = 0.8, unstranded_threshold: float = 0.1
) -> tuple[str, float, float, float]:
    """Upstream's strandedness call from salmon's ``lib_format_counts.json``."""
    def total(keys: tuple[str, ...]) -> float:
        return sum(lib_format_counts.get(k) or 0 for k in keys)

    return calculate_strandedness(
        total(("SF", "ISF", "MSF", "OSF")),
        total(("SR", "ISR", "MSR", "OSR")),
        total(("IU", "U", "MU")),
        stranded_threshold,
        unstranded_threshold,
    )


async def infer_strandedness(
    sample: Sample,
    reads_1: File,
    reads_2: File | None,
    index: Dir,
    genome: Genome,
    stranded_threshold: float = 0.8,
    unstranded_threshold: float = 0.1,
) -> StrandednessAnalysis:
    """Subsample the reads and let salmon classify the library, as upstream does."""
    sub_1, sub_2 = await fq_subsample(
        reads_1, reads_2, args="--record-count 1000000 --seed 1", prefix=f"{sample.id}.subsampled"
    )
    results = await salmon_quant_reads(
        reads_1=sub_1,
        reads_2=sub_2,
        index=index,
        gtf=genome.gtf,
        lib_type="A",
        args="--skipQuant",
    )
    counts_file = await results.get_file("lib_format_counts.json")
    if counts_file is None:
        raise RuntimeError(
            f"salmon produced no lib_format_counts.json for sample {sample.id!r} (strandedness 'auto'); "
            "check that the salmon index matches the reads, or set strandedness explicitly"
        )
    async with counts_file.open("rb") as fh:
        counts = json.loads(bytes(await fh.read()))
    inferred, fwd, rev, unstr = salmon_strandedness(counts, stranded_threshold, unstranded_threshold)
    return StrandednessAnalysis(inferred, fwd, rev, unstr, counts_file)


def salmon_lib_type(sample: Sample) -> str:
    """Map samplesheet strandedness to a salmon library type, as upstream does."""
    single = {"forward": "SF", "reverse": "SR", "unstranded": "U"}
    paired = {"forward": "ISF", "reverse": "ISR", "unstranded": "IU"}
    if sample.strandedness == "auto":
        return "A"
    return (single if sample.single_end else paired)[sample.strandedness]


async def quantify_salmon_bam(sample: Sample, alignment: StarAlignment, genome: Genome) -> Dir:
    """Salmon alignment-mode quantification of the STAR transcriptome BAM."""
    return await salmon_quant_bam(
        bam=alignment.transcriptome_bam,
        transcript_fasta=genome.transcript_fasta,
        gtf=genome.gtf,
        lib_type=salmon_lib_type(sample),
    )


@dataclass
class RunOptions:
    seq_platform: str = ""
    seq_center: str = ""
    skip_fastqc: bool = False
    skip_trimming: bool = False
    min_trimmed_reads: int = 10000
    stranded_threshold: float = 0.8
    unstranded_threshold: float = 0.1
    skip_linting: bool = False
    extra_fqlint_args: str = "--disable-validator P001"
    skip_markduplicates: bool = False
    skip_dupradar: bool = False
    skip_qualimap: bool = False
    skip_rseqc: bool = False
    rseqc_modules: tuple[str, ...] = RSEQC_MODULES
    skip_biotype_qc: bool = False
    featurecounts_group_type: str = "gene_biotype"
    featurecounts_feature_type: str = "exon"
    skip_stringtie: bool = False
    skip_bigwig: bool = False
    min_mapped_reads: float = 5.0


async def run_sample(
    sample: Sample,
    genome: "asyncio.Task[Genome]",
    strandedness_index: "asyncio.Task[Dir] | None",
    bbsplit_index_dir: "asyncio.Task[Dir] | None",
    opts: RunOptions,
) -> SampleResult | float:
    """Preprocess, (infer strandedness,) align and quantify one sample.

    Returns the surviving read count instead of a result when trimming left
    fewer than ``opts.min_trimmed_reads`` reads (the sample is dropped).
    """
    reads = await preprocess_reads(sample, opts, bbsplit_index_dir)
    if reads.reads_after_trimming is not None and reads.reads_after_trimming < opts.min_trimmed_reads:
        return reads.reads_after_trimming
    ref = await genome

    analysis = None
    if sample.strandedness == "auto":
        assert strandedness_index is not None
        analysis = await infer_strandedness(
            sample,
            reads.reads_1,
            reads.reads_2,
            await strandedness_index,
            ref,
            opts.stranded_threshold,
            opts.unstranded_threshold,
        )
        # As upstream, an undetermined library is treated as unstranded.
        inferred = "unstranded" if analysis.inferred == "undetermined" else analysis.inferred
        sample = replace(sample, strandedness=inferred)

    alignment = await align_star(sample, reads.reads_1, reads.reads_2, ref, opts.seq_platform, opts.seq_center)

    async def markdup() -> MarkedDuplicates | None:
        return None if opts.skip_markduplicates else await mark_duplicates(sample, alignment, ref)

    async def genome_bam_qc() -> tuple[
        MarkedDuplicates | None,
        DupradarResult | None,
        Dir | None,
        dict[str, Dir],
        tuple[FeatureCountsResult, BiotypeQC] | None,
        StringTieResult | None,
        dict[str, File],
    ]:
        marked = await markdup()
        # As upstream, the QC runs on the duplicate-marked BAM when there is one.
        bam = marked.bam if marked is not None else alignment.bam
        bai = marked.bai if marked is not None else alignment.bai

        async def run_dupradar() -> DupradarResult | None:
            if opts.skip_dupradar:
                return None
            return await dupradar(bam, ref.gtf, sample.id, sample.strandedness, sample.single_end)

        async def run_qualimap() -> Dir | None:
            if opts.skip_qualimap:
                return None
            # As upstream: Qualimap on a name-sorted copy, told it's sorted.
            namesorted = await samtools_sort(bam=bam, args="-n")
            return await qualimap_rnaseq(
                namesorted,
                ref.gtf,
                sample.id,
                sample.strandedness,
                sample.single_end,
                args="--sorted",
                bam_name=f"{sample.id}.namesorted.bam",  # upstream's name; MultiQC reads it from the report
            )

        async def run_rseqc() -> dict[str, Dir]:
            if opts.skip_rseqc or not opts.rseqc_modules:
                return {}
            return await rseqc(bam, bai, ref.gene_bed, sample.id, sample.single_end, opts.rseqc_modules)

        async def run_biotype_qc() -> tuple[FeatureCountsResult, BiotypeQC] | None:
            if opts.skip_biotype_qc:
                return None
            # As upstream: featureCounts grouped by biotype, then MultiQC tables.
            counts = await featurecounts(
                bam,
                ref.gtf,
                sample.id,
                sample.strandedness,
                sample.single_end,
                args=f"-B -C -g {opts.featurecounts_group_type} -t {opts.featurecounts_feature_type}",
            )
            return counts, await multiqc_custom_biotype(counts.counts, prefix=sample.id)

        async def run_stringtie() -> StringTieResult | None:
            if opts.skip_stringtie:
                return None
            # As upstream: reference-guided (-G GTF), estimating known transcripts only (-e).
            return await stringtie(bam, sample.id, gtf=ref.gtf, strandedness=sample.strandedness, args="-v -e")

        async def run_bigwig() -> dict[str, File]:
            return {} if opts.skip_bigwig else await bigwig_coverage(sample, bam, ref.chrom_sizes)

        dup, qm, rs, bt, st, bw = await asyncio.gather(
            run_dupradar(), run_qualimap(), run_rseqc(), run_biotype_qc(), run_stringtie(), run_bigwig()
        )
        return marked, dup, qm, rs, bt, st, bw

    # As upstream: samples mapping below min_mapped_reads are still quantified
    # (salmon runs on the transcriptome BAM) but skip every genome-BAM step.
    percent_mapped = await star_percent_mapped(alignment.log_final)
    passed = percent_mapped >= opts.min_mapped_reads

    async def gated_genome_bam_qc():
        return await genome_bam_qc() if passed else (None, None, None, {}, None, None, {})

    salmon, (marked, dup, qm, rs, bt, st, bw) = await asyncio.gather(
        quantify_salmon_bam(sample, alignment, ref), gated_genome_bam_qc()
    )
    return SampleResult(
        sample=sample.id,
        percent_mapped=percent_mapped,
        mapping_passed=passed,
        strandedness=sample.strandedness,
        strandedness_analysis=analysis,
        preprocessing=reads,
        alignment=alignment,
        salmon=salmon,
        markduplicates=marked,
        dupradar=dup,
        qualimap=qm,
        rseqc=rs,
        biotype_counts=bt[0] if bt else None,
        biotype_qc=bt[1] if bt else None,
        stringtie=st,
        bigwig=bw,
    )


async def write_samplesheet(samples: list[Sample]) -> File:
    """Write the samples as an upstream-style samplesheet CSV (one row per run).

    Used as the SummarizedExperiment column metadata, as upstream passes its
    input samplesheet; the ``.csv`` name tells the reader the delimiter.
    """
    out = Path(tempfile.mkdtemp(prefix="samplesheet_")) / "samplesheet.csv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        for s in samples:
            mates = s.fastq_2 or [None] * len(s.fastq_1)
            for r1, r2 in zip(s.fastq_1, mates):
                writer.writerow([s.id, r1.path, r2.path if r2 is not None else "", s.strandedness])
    return await File.from_local(out)


async def merge_quantifications(
    samples: list[Sample],
    salmon_results: dict[str, Dir],
    gtf: File,
    gene_id: str = "gene_id",
    extra_attributes: str = "gene_name",
) -> MergedQuantification:
    """Upstream's tximport + SummarizedExperiment merge of per-sample salmon results.

    tx2gene only needs one sample's quant file to discover the transcript ID
    attribute (all samples share a transcriptome), as upstream does.
    """
    first = next(iter(salmon_results.values()))
    tx2gene_tsv, quants, coldata = await asyncio.gather(
        tx2gene(gtf, first, quant_type="salmon", gene_id=gene_id, extra=extra_attributes, prefix="salmon.merged"),
        collect_quants(salmon_results),
        write_samplesheet(samples),
    )
    txi = await tximport(quants, tx2gene_tsv, quant_type="salmon", prefix="salmon.merged")
    gene_rds, transcript_rds = await asyncio.gather(
        summarized_experiment(
            {
                "counts": txi.counts_gene,
                "counts_length_scaled": txi.counts_gene_length_scaled,
                "counts_scaled": txi.counts_gene_scaled,
                "lengths": txi.lengths_gene,
                "tpm": txi.tpm_gene,
            },
            rowdata=tx2gene_tsv,
            coldata=coldata,
            prefix="salmon.merged.gene",
        ),
        summarized_experiment(
            {"counts": txi.counts_transcript, "lengths": txi.lengths_transcript, "tpm": txi.tpm_transcript},
            rowdata=tx2gene_tsv,
            coldata=coldata,
            prefix="salmon.merged.transcript",
        ),
    )
    return MergedQuantification(tx2gene=tx2gene_tsv, tximport=txi, gene_rds=gene_rds, transcript_rds=transcript_rds)


async def rnaseq(
    samples: list[Sample],
    fasta: File,
    gtf: File | None = None,
    gff: File | None = None,
    transcript_fasta: File | None = None,
    additional_fasta: File | None = None,
    salmon_index_dir: Dir | None = None,
    bbsplit_fasta_list: dict[str, File] | None = None,
    bbsplit_index_dir: Dir | None = None,
    skip_bbsplit: bool = True,
    skip_linting: bool = False,
    extra_fqlint_args: str = "--disable-validator P001",
    skip_markduplicates: bool = False,
    skip_dupradar: bool = False,
    skip_qualimap: bool = False,
    skip_rseqc: bool = False,
    rseqc_modules: tuple[str, ...] = RSEQC_MODULES,
    skip_biotype_qc: bool = False,
    featurecounts_group_type: str = "gene_biotype",
    featurecounts_feature_type: str = "exon",
    skip_stringtie: bool = False,
    skip_bigwig: bool = False,
    skip_deseq2_qc: bool = False,
    min_mapped_reads: float = 5.0,
    skip_multiqc: bool = False,
    seq_platform: str = "",
    seq_center: str = "",
    skip_fastqc: bool = False,
    skip_trimming: bool = False,
    min_trimmed_reads: int = 10000,
    stranded_threshold: float = 0.8,
    unstranded_threshold: float = 0.1,
) -> RnaseqResult:
    """Prepare the genome while preprocessing every sample, then align and quantify them in parallel.

    ``salmon_index_dir`` is only used to infer strandedness for
    ``strandedness="auto"`` samples; without it one is built from the genome
    and transcript FASTA (decoy-aware), and only if some sample is ``auto``.

    BBSplit (off by default, as upstream) removes reads that map better to
    other genomes: give ``bbsplit_fasta_list`` (name → FASTA) to build an index
    against the prepared genome, or a prebuilt ``bbsplit_index_dir``.
    """
    if not skip_bbsplit and bbsplit_fasta_list is None and bbsplit_index_dir is None:
        raise ValueError("BBSplit needs bbsplit_fasta_list or bbsplit_index_dir (or skip_bbsplit=True)")
    ids = [s.id for s in samples]
    if len(ids) != len(set(ids)):
        raise ValueError("sample ids must be unique; put a sample's runs in one Sample")

    genome = asyncio.create_task(
        prepare_genome(
            fasta=fasta,
            gtf=gtf,
            gff=gff,
            transcript_fasta=transcript_fasta,
            additional_fasta=additional_fasta,
        )
    )
    opts = RunOptions(
        seq_platform=seq_platform,
        seq_center=seq_center,
        skip_fastqc=skip_fastqc,
        skip_trimming=skip_trimming,
        min_trimmed_reads=min_trimmed_reads,
        stranded_threshold=stranded_threshold,
        unstranded_threshold=unstranded_threshold,
        skip_linting=skip_linting,
        extra_fqlint_args=extra_fqlint_args,
        skip_markduplicates=skip_markduplicates,
        skip_dupradar=skip_dupradar,
        skip_qualimap=skip_qualimap,
        skip_rseqc=skip_rseqc,
        rseqc_modules=tuple(rseqc_modules),
        skip_biotype_qc=skip_biotype_qc,
        featurecounts_group_type=featurecounts_group_type,
        featurecounts_feature_type=featurecounts_feature_type,
        skip_stringtie=skip_stringtie,
        skip_bigwig=skip_bigwig,
        min_mapped_reads=min_mapped_reads,
    )

    async def strandedness_index() -> Dir:
        if salmon_index_dir is not None:
            return salmon_index_dir
        ref = await genome
        return await salmon_index(transcript_fasta=ref.transcript_fasta, genome_fasta=ref.fasta)

    index = asyncio.create_task(strandedness_index()) if any(s.strandedness == "auto" for s in samples) else None

    async def build_bbsplit_index() -> Dir:
        if bbsplit_index_dir is not None:
            return bbsplit_index_dir
        assert bbsplit_fasta_list is not None
        ref = await genome
        return await bbsplit_index(ref.fasta, bbsplit_fasta_list, args=BBSPLIT_ARGS)

    split_index = None if skip_bbsplit else asyncio.create_task(build_bbsplit_index())
    try:
        outcomes = await asyncio.gather(*(run_sample(s, genome, index, split_index, opts) for s in samples))
        ref = await genome
    finally:
        for task in (genome, index, split_index):
            if task is not None and not task.done():
                task.cancel()
    passed = [o for o in outcomes if isinstance(o, SampleResult)]
    failed = {s.id: o for s, o in zip(samples, outcomes) if not isinstance(o, SampleResult)}
    if not passed:
        raise RuntimeError(f"no samples passed the {min_trimmed_reads}-read trimming threshold: {failed}")
    merged = await merge_quantifications(samples, {r.sample: r.salmon for r in passed}, ref.gtf)
    # As upstream: DESeq2 QC of the merged length-scaled gene counts, labelled by aligner.
    dqc = None if skip_deseq2_qc else await deseq2_qc(merged.tximport.counts_gene_length_scaled, label="star_salmon")
    result = RnaseqResult(genome=ref, samples=passed, salmon=merged, failed_trimming=failed, deseq2_qc=dqc)
    if not skip_multiqc:
        result.multiqc = await multiqc_report(
            samples,
            result,
            min_trimmed_reads=min_trimmed_reads,
            rseqc_modules=() if skip_rseqc else tuple(rseqc_modules),
            stranded_threshold=stranded_threshold,
            unstranded_threshold=unstranded_threshold,
        )
    return result


# ---------------------------------------------------------------------------
# MultiQC
# ---------------------------------------------------------------------------


def round_one_decimal(v: float | None) -> float | None:
    """Groovy's Math.round(v * 10) / 10 (half-up), as upstream rounds report values."""
    return None if v is None else math.floor(v * 10 + 0.5) / 10.0


def infer_experiment_strandedness(text: str, stranded_threshold: float, unstranded_threshold: float):
    """Upstream's getInferexperimentStrandedness: RSeQC fractions -> (strandedness, fwd %, rev %, unstr %)."""
    fwd = rev = unstr = 0.0
    patterns = {
        "unstr": r"Fraction of reads failed to determine:\s([\d\.]+)",
        "fwd_se": r'Fraction of reads explained by "\++,--":\s([\d\.]+)',
        "rev_se": r'Fraction of reads explained by "\+-,-\+":\s([\d\.]+)',
        "fwd_pe": r'Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)',
        "rev_pe": r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)',
    }
    for line in text.splitlines():
        for key, pattern in patterns.items():
            if m := re.search(pattern, line):
                value = float(m.group(1)) * 100
                if key == "unstr":
                    unstr = value
                elif key.startswith("fwd"):
                    fwd = value
                else:
                    rev = value
    return calculate_strandedness(fwd, rev, unstr, stranded_threshold, unstranded_threshold)


@dataclass
class StrandCall:
    inferred: str
    forward: float
    reverse: float
    unstranded: float

    def certainty(self) -> float | None:
        stranded = self.forward + self.reverse
        if stranded == 0:
            return None
        if self.inferred == "forward":
            return self.forward / stranded * 100
        if self.inferred == "reverse":
            return self.reverse / stranded * 100
        return None


def load_multiqc_asset(name: str) -> dict:
    """Upstream's loadMultiqcAsset: the YAML (merge keys resolved) without `_`-prefixed keys."""
    parsed = yaml.safe_load(asset(name).read_text())
    return {k: v for k, v in parsed.items() if not str(k).startswith("_")}


def strand_check_summary(rows: list[tuple[str, str, str, StrandCall | None, StrandCall | None]]) -> str:
    """Upstream's strandCheckSummaryYaml, as MultiQC custom-content JSON."""
    config = load_multiqc_asset("strand_check_summary.yaml")
    header_keys = list(config["headers"])
    data = {}
    for sample, provided, status, salmon, rseqc_call in sorted(rows, key=lambda r: r[0]):
        raw = {
            "provided": provided,
            "salmon_inferred": salmon.inferred if salmon else "-",
            "salmon_pct": round_one_decimal(salmon.certainty()) if salmon else None,
            "salmon_s": round_one_decimal(salmon.forward) if salmon else None,
            "salmon_a": round_one_decimal(salmon.reverse) if salmon else None,
            "salmon_u": round_one_decimal(salmon.unstranded) if salmon else None,
            "rseqc_inferred": rseqc_call.inferred if rseqc_call else "-",
            "rseqc_pct": round_one_decimal(rseqc_call.certainty()) if rseqc_call else None,
            "rseqc_s": round_one_decimal(rseqc_call.forward) if rseqc_call else None,
            "rseqc_a": round_one_decimal(rseqc_call.reverse) if rseqc_call else None,
            "rseqc_u": round_one_decimal(rseqc_call.unstranded) if rseqc_call else None,
            "status": status,
        }
        unknown = set(raw) - set(header_keys)
        if unknown:
            raise ValueError(f"strand_check_summary.yaml headers do not declare columns: {unknown}")
        data[sample] = {k: raw[k] for k in header_keys if raw[k] is not None}
    return json.dumps({**config, "data": data}, indent=4)


def strand_check_composition(rows: list[tuple[str, str, str, StrandCall | None, StrandCall | None]]) -> str:
    """Upstream's strandCheckCompositionYaml, as MultiQC custom-content JSON."""
    summary = load_multiqc_asset("strand_check_summary.yaml")
    config = {
        **load_multiqc_asset("strand_check_composition.yaml"),
        **{k: summary[k] for k in ("parent_id", "parent_name", "parent_description")},
    }

    def composition(call: StrandCall) -> dict:
        return {
            "Sense": round_one_decimal(call.forward),
            "Antisense": round_one_decimal(call.reverse),
            "Unstranded": round_one_decimal(call.unstranded),
        }

    rseqc_data, salmon_data = {}, {}
    for sample, _provided, _status, salmon, rseqc_call in sorted(rows, key=lambda r: r[0]):
        if rseqc_call:
            rseqc_data[sample] = composition(rseqc_call)
        if salmon:
            salmon_data[sample] = composition(salmon)
    datasets, labels = [], []
    if rseqc_data:
        datasets.append(rseqc_data)
        labels.append("RSeQC")
    if salmon_data:
        datasets.append(salmon_data)
        labels.append("Salmon")
    pconfig = dict(config["pconfig"])
    if len(datasets) > 1:
        pconfig["data_labels"] = [{"name": label, "ylab": pconfig["ylab"]} for label in labels]
    config["pconfig"] = pconfig
    config["data"] = datasets[0] if len(datasets) == 1 else datasets
    return json.dumps(config, indent=4)


def sample_merge_yaml(samples: list[Sample]) -> str:
    """Upstream's multiqcSampleMergeYaml: merge `<sample>_1` / `<sample>_2` rows of paired-end samples."""
    ids = sorted({s.id for s in samples if not s.single_end})
    if not ids:
        return "table_sample_merge: {}\n"

    def pattern(sample_id: str, read: int) -> str:
        esc = re.sub(r"[\\^$.|?*+()\[\]{}/]", lambda m: "\\" + m.group(0), sample_id).replace("'", "''")
        return f"    - type: regex\n      pattern: '(?<=^{esc})_{read}$'"

    r1 = "\n".join(pattern(i, 1) for i in ids)
    r2 = "\n".join(pattern(i, 2) for i in ids)
    return f'table_sample_merge:\n  "Read 1":\n{r1}\n  "Read 2":\n{r2}\n'


def name_replacements(samples: list[Sample]) -> str:
    """Upstream's multiqcNameReplacements: FASTQ simple names -> sample ids, where they differ."""

    def simple_name(f: File) -> str:
        return f.path.rstrip("/").rsplit("/", 1)[-1].split(".", 1)[0]

    lines = []
    for s in samples:
        suffixes = ("", "") if s.single_end else ("_1", "_2")
        first = simple_name(s.fastq_1[0])
        if first != s.id:
            lines.append(f"{first}\t{s.id}{suffixes[0]}")
            if not s.single_end:
                lines.append(f"{simple_name(s.fastq_2[0])}\t{s.id}{suffixes[1]}")
    return "".join(f"{line}\n" for line in lines)


async def multiqc_report(
    samples: list[Sample],
    result: "RnaseqResult",
    min_trimmed_reads: int,
    rseqc_modules: tuple[str, ...],
    stranded_threshold: float,
    unstranded_threshold: float,
) -> MultiQCResult:
    """Lay every QC output out under upstream's file names, add upstream's custom content, run MultiQC."""
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

    def write(rel: str, text: str) -> None:
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(text)

    strand_rows = []
    for r in result.samples:
        sid, pre, a = r.sample, r.preprocessing, r.alignment
        await put_dir(pre.raw_fastqc, f"fastqc/raw/{sid}", ("_fastqc.zip",))
        if pre.trimming is not None:
            for report in pre.trimming.reports:
                put(report, f"trimgalore/{sid}/{report.path.rsplit('/', 1)[-1]}")
            await put_dir(pre.trimming.results, f"fastqc/trim/{sid}", ("_fastqc.zip",))
        if pre.bbsplit is not None:
            put(pre.bbsplit.stats, f"bbsplit/{sid}.stats.txt")
        put(a.log_final, f"star/{sid}.Log.final.out")
        for kind, f in (("stats", a.stats), ("flagstat", a.flagstat), ("idxstats", a.idxstats)):
            put(f, f"samtools/{sid}.sorted.bam.{kind}")
        md = r.markduplicates
        if md is not None:
            put(md.metrics, f"picard/{sid}.markdup.sorted.metrics.txt")
            for kind, f in (("stats", md.stats), ("flagstat", md.flagstat), ("idxstats", md.idxstats)):
                put(f, f"samtools/{sid}.markdup.sorted.bam.{kind}")
        await put_dir(r.salmon, f"salmon/{sid}")
        await put_dir(r.qualimap, f"qualimap/{sid}")
        if r.dupradar is not None:
            put(r.dupradar.intercept_mqc, f"dupradar/{sid}_dup_intercept_mqc.txt")
            put(r.dupradar.curve_mqc, f"dupradar/{sid}_duprateExpDensCurve_mqc.txt")
        for module, d in r.rseqc.items():
            await put_dir(d, f"rseqc/{module}", (".txt", ".log", ".r", ".xls"))
        if r.biotype_qc is not None:
            put(r.biotype_qc.counts_mqc, f"biotype/{sid}.biotype_counts_mqc.tsv")
            put(r.biotype_qc.rrna_mqc, f"biotype/{sid}.biotype_counts_rrna_mqc.tsv")

        # Strandedness checks: classify against RSeQC infer_experiment when it ran, else
        # surface Salmon's call for 'auto' samples (upstream's classifyStrand logic).
        a_ = r.strandedness_analysis
        salmon = StrandCall(a_.inferred, a_.forward_pct, a_.reverse_pct, a_.unstranded_pct) if a_ else None
        infer = r.rseqc.get("infer_experiment")
        if "infer_experiment" in rseqc_modules and infer is not None:
            txt = await infer.get_file(f"{sid}.infer_experiment.txt")
            assert txt is not None
            async with txt.open("rb") as fh:
                call = StrandCall(
                    *infer_experiment_strandedness(
                        bytes(await fh.read()).decode(), stranded_threshold, unstranded_threshold
                    )
                )
            if salmon is not None:
                provided = "auto"
                status = "pass" if salmon.inferred == call.inferred and call.inferred != "undetermined" else "fail"
            else:
                provided = r.strandedness
                status = "pass" if r.strandedness == call.inferred else "fail"
            strand_rows.append((sid, provided, status, salmon, call))
        elif salmon is not None:
            strand_rows.append((sid, "auto", "-", salmon, None))

    if result.deseq2_qc is not None:
        put(result.deseq2_qc.pca_multiqc, "deseq2/star_salmon.pca.vals_mqc.tsv")
        put(result.deseq2_qc.dists_multiqc, "deseq2/star_salmon.sample.dists_mqc.tsv")

    # Upstream lists samples at or below min_trimmed_reads (the filter itself drops < only).
    trimmed = dict(result.failed_trimming)
    for r in result.samples:
        n = r.preprocessing.reads_after_trimming
        if n is not None and n <= min_trimmed_reads:
            trimmed[r.sample] = n
    if trimmed:
        rows = "".join(f"{sid}\t{n}\n" for sid, n in trimmed.items())
        write("custom/fail_trimmed_samples_mqc.tsv", "Sample\tReads after trimming\n" + rows)
    failed_mapping = [r for r in result.samples if not r.mapping_passed]
    if failed_mapping:
        rows = "".join(f"{r.sample}\t{r.percent_mapped}\n" for r in failed_mapping)
        write(
            "custom/fail_mapped_samples_mqc.tsv",
            asset("sample_status_header.txt").read_text() + "Sample\tSTAR uniquely mapped reads (%)\n" + rows,
        )
    if strand_rows:
        write("custom/strand_check_summary_mqc.json", strand_check_summary(strand_rows))
        write("custom/strand_check_composition_mqc.json", strand_check_composition(strand_rows))

    await asyncio.gather(*downloads)
    data = await Dir.from_local(root)

    config_dir = Path(tempfile.mkdtemp(prefix="multiqc_cfg_"))
    (config_dir / "multiqc_sample_merge.yml").write_text(sample_merge_yaml(samples))
    configs = [
        await File.from_local(str(asset("multiqc_config.yml"))),
        await File.from_local(config_dir / "multiqc_sample_merge.yml"),
    ]
    replacements = name_replacements(samples)
    replace_file = None
    if replacements:
        (config_dir / "name_replacement.txt").write_text(replacements)
        replace_file = await File.from_local(config_dir / "name_replacement.txt")
    return await multiqc(data, configs=configs, replace_names=replace_file, prefix="multiqc_report")
