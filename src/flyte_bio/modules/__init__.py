"""flyte_bio.modules — typed shell-task wrappers around bio CLI tools.

Each submodule wraps one tool family, sharing a single biocontainer image
across its tasks:

- :mod:`flyte_bio.modules.bbmap` — BBSplit (index build + read binning).
- :mod:`flyte_bio.modules.bedtools` — genome arithmetic (intersect, sort,
  merge).
- :mod:`flyte_bio.modules.cat` — file concatenation (``cat_fastq``).
- :mod:`flyte_bio.modules.dupradar` — duplication rate vs expression QC.
- :mod:`flyte_bio.modules.fastqc` — read quality-control reports.
- :mod:`flyte_bio.modules.catadditionalfasta` — append an extra FASTA + GTF.
- :mod:`flyte_bio.modules.fq` — FASTQ subsampling and linting.
- :mod:`flyte_bio.modules.gffread` — GFF/GTF conversion + transcript FASTA.
- :mod:`flyte_bio.modules.gtf2bed` — derive a BED12 gene model from a GTF.
- :mod:`flyte_bio.modules.gtffilter` — restrict a GTF to a genome's sequences.
- :mod:`flyte_bio.modules.gunzip` — single-file gzip decompression.
- :mod:`flyte_bio.modules.picard` — Picard MarkDuplicates.
- :mod:`flyte_bio.modules.qualimap` — RNA-seq alignment QC.
- :mod:`flyte_bio.modules.rseqc` — RSeQC BAM QC (7 scripts).
- :mod:`flyte_bio.modules.salmon` — transcript quantification (index, quant).
- :mod:`flyte_bio.modules.samtools` — alignment sort/index/stats/faidx.
- :mod:`flyte_bio.modules.star` — spliced RNA-seq aligner (index, align).
- :mod:`flyte_bio.modules.summarizedexperiment` — bundle matrices into an RDS.
- :mod:`flyte_bio.modules.trimgalore` — adapter/quality trimming.
- :mod:`flyte_bio.modules.tx2gene` — transcript → gene table from a GTF.
- :mod:`flyte_bio.modules.tximport` — count/TPM matrices from quantifications.
- :mod:`flyte_bio.modules.untar` — tar archive extraction.

The module-level :data:`env` here is an aggregate
:class:`flyte.TaskEnvironment` depending on every submodule's env.
Pipelines depend on it once to gain access to every wrapped tool::

    from flyte_bio.modules import env as modules_env
"""



import flyte

from .bbmap import env as bbmap_env
from .bedtools import env as bedtools_env
from .cat import env as cat_env
from .catadditionalfasta import env as catadditionalfasta_env
from .dupradar import env as dupradar_env
from .fastqc import env as fastqc_env
from .fq import env as fq_env
from .gffread import env as gffread_env
from .gtf2bed import env as gtf2bed_env
from .gtffilter import env as gtffilter_env
from .gunzip import env as gunzip_env
from .picard import env as picard_env
from .qualimap import env as qualimap_env
from .rseqc import env as rseqc_env
from .salmon import env as salmon_env
from .samtools import env as samtools_env
from .star import env as star_env
from .summarizedexperiment import env as summarizedexperiment_env
from .trimgalore import env as trimgalore_env
from .tx2gene import env as tx2gene_env
from .tximport import env as tximport_env
from .untar import env as untar_env

env = flyte.TaskEnvironment(
    name="flyte_bio_modules",
    depends_on=[
        bbmap_env,
        bedtools_env,
        cat_env,
        catadditionalfasta_env,
        dupradar_env,
        fastqc_env,
        fq_env,
        gffread_env,
        gtf2bed_env,
        gtffilter_env,
        gunzip_env,
        picard_env,
        qualimap_env,
        rseqc_env,
        salmon_env,
        samtools_env,
        star_env,
        summarizedexperiment_env,
        trimgalore_env,
        tx2gene_env,
        tximport_env,
        untar_env,
    ],
)

__all__ = [
    "bbmap_env",
    "bedtools_env",
    "cat_env",
    "catadditionalfasta_env",
    "dupradar_env",
    "env",
    "fastqc_env",
    "fq_env",
    "gffread_env",
    "gtf2bed_env",
    "gtffilter_env",
    "gunzip_env",
    "picard_env",
    "qualimap_env",
    "rseqc_env",
    "salmon_env",
    "samtools_env",
    "star_env",
    "summarizedexperiment_env",
    "trimgalore_env",
    "tx2gene_env",
    "tximport_env",
    "untar_env",
]
