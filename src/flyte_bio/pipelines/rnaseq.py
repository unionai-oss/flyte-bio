"""rnaseq — RNA-seq quantification (STAR alignment + salmon), rnaseq 3.26.0.

Composes :mod:`flyte_bio.modules` tasks into the ``star_salmon`` path:

- :func:`prepare_genome` — decompress references, filter the GTF to the
  genome's sequences, append an additional FASTA, derive the gene BED,
  transcript FASTA and chrom sizes, and build the STAR index.
- :func:`align_star` — merge a sample's runs, align with STAR (emitting a
  transcriptome BAM), then sort / index / stats the genome BAM.
- :func:`quantify_salmon_bam` — salmon alignment-mode quantification of the
  transcriptome BAM.
- :func:`rnaseq` — run all of the above, samples in parallel.

These are plain async functions, not tasks: they run inside the caller's
task and fan out to the module tasks, so the caller's
:class:`flyte.TaskEnvironment` must ``depends_on`` :data:`flyte_bio.modules.env`.

Not yet covered (the rest of the upstream default path): read QC/trimming,
bbsplit, rRNA removal, strandedness inference, UMI deduplication, duplicate
marking, tximport/SummarizedExperiment merging and the QC reports. Samples
with ``strandedness="auto"`` are quantified with salmon's own library-type
auto-detection (``A``) instead of upstream's subsampled inference.
"""

import asyncio
from dataclasses import dataclass, field

from flyte.io import Dir, File

from flyte_bio.modules.cat import cat_fastq
from flyte_bio.modules.catadditionalfasta import cat_additional_fasta
from flyte_bio.modules.gffread import gffread_gff_to_gtf, gffread_transcripts_fasta
from flyte_bio.modules.gtf2bed import gtf2bed
from flyte_bio.modules.gtffilter import gtf_filter
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.salmon import salmon_quant_bam
from flyte_bio.modules.samtools import (
    samtools_faidx,
    samtools_flagstat,
    samtools_idxstats,
    samtools_index,
    samtools_sort,
    samtools_stats,
)
from flyte_bio.modules.star import star_align, star_genome_generate

STRANDEDNESS = ("auto", "forward", "reverse", "unstranded")


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
class SampleResult:
    sample: str
    alignment: StarAlignment
    salmon: Dir


@dataclass
class RnaseqResult:
    genome: Genome
    samples: list[SampleResult]


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


async def align_star(sample: Sample, genome: Genome, seq_platform: str = "", seq_center: str = "") -> StarAlignment:
    """Merge the sample's runs, align with STAR, then sort / index / stats the BAM."""
    if any(c.isspace() for c in sample.id):
        raise ValueError(f"sample id {sample.id!r} must not contain whitespace")

    reads_1, reads_2 = await asyncio.gather(
        cat_fastq(sample.fastq_1),
        cat_fastq(sample.fastq_2) if sample.fastq_2 else asyncio.sleep(0, result=None),
    )
    output = await star_align(
        reads_1=reads_1,
        reads_2=[reads_2] if reads_2 is not None else [],
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


async def run_sample(sample: Sample, genome: Genome, seq_platform: str, seq_center: str) -> SampleResult:
    alignment = await align_star(sample, genome, seq_platform, seq_center)
    salmon = await quantify_salmon_bam(sample, alignment, genome)
    return SampleResult(sample=sample.id, alignment=alignment, salmon=salmon)


async def rnaseq(
    samples: list[Sample],
    fasta: File,
    gtf: File | None = None,
    gff: File | None = None,
    transcript_fasta: File | None = None,
    additional_fasta: File | None = None,
    seq_platform: str = "",
    seq_center: str = "",
) -> RnaseqResult:
    """Prepare the genome once, then align and quantify every sample in parallel."""
    ids = [s.id for s in samples]
    if len(ids) != len(set(ids)):
        raise ValueError("sample ids must be unique; put a sample's runs in one Sample")

    genome = await prepare_genome(
        fasta=fasta,
        gtf=gtf,
        gff=gff,
        transcript_fasta=transcript_fasta,
        additional_fasta=additional_fasta,
    )
    results = await asyncio.gather(*(run_sample(s, genome, seq_platform, seq_center) for s in samples))
    return RnaseqResult(genome=genome, samples=list(results))
