"""flyte_bio.modules — typed shell-task wrappers around bio CLI tools.

Each submodule wraps one tool family, sharing a single biocontainer image
across its tasks:

- :mod:`flyte_bio.modules.amulety` — AIRR translation and sequence embeddings.
- :mod:`flyte_bio.modules.anndata` — subset an h5ad to listed barcodes.
- :mod:`flyte_bio.modules.anndatar` — h5ad to Seurat / SingleCellExperiment RDS.
- :mod:`flyte_bio.modules.aria2` — file downloads.
- :mod:`flyte_bio.modules.art` — Illumina read simulation.
- :mod:`flyte_bio.modules.bakta` — bacterial genome annotation (database download, annotate).
- :mod:`flyte_bio.modules.bbmap` — BBSplit (index build + read binning).
- :mod:`flyte_bio.modules.bcftools` — VCF/BCF utilities (stats, sort, norm, view, filter, annotate, merge, mpileup).
- :mod:`flyte_bio.modules.bedops` — GTF to BED conversion.
- :mod:`flyte_bio.modules.bedtools` — genome arithmetic (intersect, sort,
  merge, genomecov).
- :mod:`flyte_bio.modules.blast` — BLAST databases and blastn search.
- :mod:`flyte_bio.modules.bowtie` — Bowtie 1 short-read alignment (build, align).
- :mod:`flyte_bio.modules.bowtie2` — Bowtie 2 short-read alignment (build, align).
- :mod:`flyte_bio.modules.bracken` — abundance re-estimation from Kraken2 reports.
- :mod:`flyte_bio.modules.busco` — assembly and proteome completeness (lineage download, assess).
- :mod:`flyte_bio.modules.bwa` — BWA short-read alignment (index).
- :mod:`flyte_bio.modules.bwamem2` — BWA-MEM2 short-read alignment (index, mem).
- :mod:`flyte_bio.modules.cellranger` — 10x Cell Ranger V(D)J (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.cat` — file concatenation (``cat_fastq``, ``cat_cat``, ``concatenate``).
- :mod:`flyte_bio.modules.deepcell` — Mesmer cell segmentation (non-commercial; in ``licensed_env``).
- :mod:`flyte_bio.modules.deseq2_qc` — DESeq2 PCA / sample-distance QC.
- :mod:`flyte_bio.modules.dupradar` — duplication rate vs expression QC.
- :mod:`flyte_bio.modules.fastp` — read trimming with QC reports.
- :mod:`flyte_bio.modules.fastqc` — read quality-control reports.
- :mod:`flyte_bio.modules.catadditionalfasta` — append an extra FASTA + GTF.
- :mod:`flyte_bio.modules.cellbender` — ambient RNA and empty droplet removal (CPU or GPU).
- :mod:`flyte_bio.modules.cellpose` — cell segmentation of microscopy images.
- :mod:`flyte_bio.modules.chromsizes` — chromosome sizes and FASTA index.
- :mod:`flyte_bio.modules.clustalo` — multiple sequence alignment and guide trees.
- :mod:`flyte_bio.modules.cooler` — Hi-C contact matrices.
- :mod:`flyte_bio.modules.csvtk` — table joins and concatenation.
- :mod:`flyte_bio.modules.cutadapt` — adapter, primer and quality trimming.
- :mod:`flyte_bio.modules.damageprofiler` — ancient-DNA damage patterns.
- :mod:`flyte_bio.modules.deeptools` — coverage matrices, heatmaps and fingerprint plots.
- :mod:`flyte_bio.modules.deepvariant` — germline small-variant calling.
- :mod:`flyte_bio.modules.diamond` — protein sequence clustering.
- :mod:`flyte_bio.modules.dragmap` — DRAGMAP short-read alignment (hash table, align).
- :mod:`flyte_bio.modules.duckdb` — CSV/TSV to Parquet conversion.
- :mod:`flyte_bio.modules.ensemblvep` — variant annotation with Ensembl VEP (cache download, annotate).
- :mod:`flyte_bio.modules.evigene` — transcript set reduction (tr2aacds).
- :mod:`flyte_bio.modules.fastqdl` — ENA/SRA FASTQ downloads.
- :mod:`flyte_bio.modules.fastqutils` — FASTQ validation.
- :mod:`flyte_bio.modules.filtlong` — long-read quality filtering.
- :mod:`flyte_bio.modules.flye` — long-read de novo assembly.
- :mod:`flyte_bio.modules.foldseek` — protein structure database building and search.
- :mod:`flyte_bio.modules.fq` — FASTQ subsampling and linting.
- :mod:`flyte_bio.modules.gatk4` — GATK tools (preprocessing, germline and somatic calling, filtering).
- :mod:`flyte_bio.modules.gawk` — run an awk program over files.
- :mod:`flyte_bio.modules.gffread` — GFF/GTF conversion + transcript FASTA.
- :mod:`flyte_bio.modules.gtf2bed` — derive a BED12 gene model from a GTF.
- :mod:`flyte_bio.modules.gtffilter` — restrict a GTF to a genome's sequences.
- :mod:`flyte_bio.modules.gtfgenefilter` — keep GTF lines on a genome's sequences (scRNA-seq).
- :mod:`flyte_bio.modules.gunzip` — single-file gzip decompression.
- :mod:`flyte_bio.modules.hisat2` — splice-aware short-read alignment (splice sites, build, align).
- :mod:`flyte_bio.modules.hmmer` — profile HMM search (hmmsearch).
- :mod:`flyte_bio.modules.htslib` — bgzip + tabix.
- :mod:`flyte_bio.modules.ilastik` — pixel classification and multicut segmentation.
- :mod:`flyte_bio.modules.interproscan` — protein domain and family annotation.
- :mod:`flyte_bio.modules.intervals` — scatter-gather genomic intervals.
- :mod:`flyte_bio.modules.kallisto` — kallisto index build.
- :mod:`flyte_bio.modules.kraken2` — taxonomic classification of reads.
- :mod:`flyte_bio.modules.krona` — interactive taxonomy charts.
- :mod:`flyte_bio.modules.macs2` — ChIP/ATAC-seq peak calling.
- :mod:`flyte_bio.modules.mafft` — multiple sequence alignment.
- :mod:`flyte_bio.modules.md5sum` — MD5 checksums of files.
- :mod:`flyte_bio.modules.megahit` — short-read and metagenome assembly.
- :mod:`flyte_bio.modules.mindagap` — tile-gap filling and duplicate spot finding.
- :mod:`flyte_bio.modules.minimap2` — long- and short-read alignment (index, align).
- :mod:`flyte_bio.modules.mmseqs` — MMseqs2 databases, clustering, taxonomy and tables.
- :mod:`flyte_bio.modules.mosdepth` — read-depth summaries.
- :mod:`flyte_bio.modules.msisensorpro` — microsatellite site scanning.
- :mod:`flyte_bio.modules.multiqc` — aggregate QC report.
- :mod:`flyte_bio.modules.multiqccustombiotype` — biotype counts for MultiQC.
- :mod:`flyte_bio.modules.nacho` — NanoString QC and normalisation.
- :mod:`flyte_bio.modules.nanolyse` — lambda phage read removal.
- :mod:`flyte_bio.modules.nanoplot` — long-read QC plots and statistics.
- :mod:`flyte_bio.modules.ncbigenomedownload` — NCBI genome downloads.
- :mod:`flyte_bio.modules.optitype` — HLA class I genotyping.
- :mod:`flyte_bio.modules.pcgr` — PCGR reference bundle download.
- :mod:`flyte_bio.modules.picard` — Picard MarkDuplicates.
- :mod:`flyte_bio.modules.pigz` — parallel gzip compression and decompression.
- :mod:`flyte_bio.modules.porechop` — Nanopore adapter trimming.
- :mod:`flyte_bio.modules.preseq` — library complexity extrapolation.
- :mod:`flyte_bio.modules.prodigal` — prokaryotic gene prediction.
- :mod:`flyte_bio.modules.prokka` — prokaryotic and viral genome annotation.
- :mod:`flyte_bio.modules.pydamage` — per-contig ancient-DNA damage estimates.
- :mod:`flyte_bio.modules.qcatch` — cell-level QC of a simpleaf quantification.
- :mod:`flyte_bio.modules.qualimap` — RNA-seq alignment QC.
- :mod:`flyte_bio.modules.quarto` — render parameterised notebooks.
- :mod:`flyte_bio.modules.quast` — assembly quality assessment.
- :mod:`flyte_bio.modules.rclone` — copy and verify data between storage remotes.
- :mod:`flyte_bio.modules.rsem` — RSEM reference preparation.
- :mod:`flyte_bio.modules.rseqc` — RSeQC BAM QC (7 scripts).
- :mod:`flyte_bio.modules.rtgtools` — RTG Tools SDF conversion.
- :mod:`flyte_bio.modules.s4pred` — protein secondary structure prediction.
- :mod:`flyte_bio.modules.salmon` — transcript quantification (index, quant).
- :mod:`flyte_bio.modules.sam2lca` — lowest-common-ancestor taxonomic assignment.
- :mod:`flyte_bio.modules.scanpy` — AnnData conversions (simpleaf to h5ad, concatenation).
- :mod:`flyte_bio.modules.simpleaf` — single-cell quantification with alevin-fry (index, quant).
- :mod:`flyte_bio.modules.samtools` — alignment sort/index/stats/faidx.
- :mod:`flyte_bio.modules.seqfu` — sequence-length statistics.
- :mod:`flyte_bio.modules.seqkit` — FASTA/FASTQ statistics and pattern replacement.
- :mod:`flyte_bio.modules.seqtk` — FASTQ subsampling.
- :mod:`flyte_bio.modules.snpeff` — variant annotation with snpEff (database download, annotate).
- :mod:`flyte_bio.modules.sortmerna` — rRNA read filtering.
- :mod:`flyte_bio.modules.sourcepredict` — metagenome source prediction.
- :mod:`flyte_bio.modules.sourmash` — MinHash sketches.
- :mod:`flyte_bio.modules.spades` — de novo genome, RNA, viral and metagenome assembly.
- :mod:`flyte_bio.modules.sratools` — SRA downloads and FASTQ dumps.
- :mod:`flyte_bio.modules.spaceranger` — 10x Space Ranger (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.star` — spliced RNA-seq aligner (index, align).
- :mod:`flyte_bio.modules.stardist` — nucleus segmentation.
- :mod:`flyte_bio.modules.stringtie` — transcript assembly / quantification.
- :mod:`flyte_bio.modules.subread` — featureCounts.
- :mod:`flyte_bio.modules.summarizedexperiment` — bundle matrices into an RDS.
- :mod:`flyte_bio.modules.thermorawfileparser` — Thermo RAW conversion.
- :mod:`flyte_bio.modules.toulligqc` — Nanopore run QC reports.
- :mod:`flyte_bio.modules.trimgalore` — adapter/quality trimming.
- :mod:`flyte_bio.modules.trimmomatic` — Illumina adapter and quality trimming.
- :mod:`flyte_bio.modules.trinity` — de novo transcriptome assembly.
- :mod:`flyte_bio.modules.trust4` — TCR/BCR repertoire reconstruction.
- :mod:`flyte_bio.modules.tx2gene` — transcript → gene table from a GTF.
- :mod:`flyte_bio.modules.tximport` — count/TPM matrices from quantifications.
- :mod:`flyte_bio.modules.ucsc` — bedClip / bedGraphToBigWig.
- :mod:`flyte_bio.modules.umicollapse` — UMI deduplication.
- :mod:`flyte_bio.modules.umitools` — UMI extraction and deduplication.
- :mod:`flyte_bio.modules.untar` — tar archive extraction.
- :mod:`flyte_bio.modules.unzip` — zip archive extraction.
- :mod:`flyte_bio.modules.vcftools` — VCF summaries (Ts/Tv, filters).
- :mod:`flyte_bio.modules.vsearch` — sequence clustering.
- :mod:`flyte_bio.modules.xz` — xz decompression.
- :mod:`flyte_bio.modules.yara` — read mapping (index, mapper).

The module-level :data:`env` here is an aggregate
:class:`flyte.TaskEnvironment` depending on every submodule's env.
Pipelines depend on it once to gain access to every wrapped tool::

    from flyte_bio.modules import env as modules_env
"""

import flyte

from .amulety import env as amulety_env
from .anndata import env as anndata_env
from .anndatar import env as anndatar_env
from .aria2 import env as aria2_env
from .art import env as art_env
from .bakta import env as bakta_env
from .bbmap import env as bbmap_env
from .bcftools import env as bcftools_env
from .bedops import env as bedops_env
from .bedtools import env as bedtools_env
from .bedtools import genomecov_env
from .blast import env as blast_env
from .bowtie import env as bowtie_env
from .bowtie2 import env as bowtie2_env
from .bracken import env as bracken_env
from .busco import env as busco_env
from .bwa import env as bwa_env
from .bwamem2 import env as bwamem2_env
from .cat import env as cat_env
from .catadditionalfasta import env as catadditionalfasta_env
from .cellbender import env as cellbender_env
from .cellpose import env as cellpose_env
from .cellranger import env as cellranger_env
from .chromsizes import env as chromsizes_env
from .clustalo import env as clustalo_env
from .cooler import env as cooler_env
from .csvtk import env as csvtk_env
from .cutadapt import env as cutadapt_env
from .damageprofiler import env as damageprofiler_env
from .deepcell import env as deepcell_env
from .deeptools import env as deeptools_env
from .deepvariant import env as deepvariant_env
from .deseq2_qc import env as deseq2_qc_env
from .diamond import env as diamond_env
from .dragmap import env as dragmap_env
from .duckdb import env as duckdb_env
from .dupradar import env as dupradar_env
from .ensemblvep import env as ensemblvep_env
from .evigene import env as evigene_env
from .fastp import env as fastp_env
from .fastqc import env as fastqc_env
from .fastqdl import env as fastqdl_env
from .fastqutils import env as fastqutils_env
from .filtlong import env as filtlong_env
from .flye import env as flye_env
from .foldseek import env as foldseek_env
from .fq import env as fq_env
from .gatk4 import env as gatk4_env
from .gawk import env as gawk_env
from .gffread import env as gffread_env
from .gtf2bed import env as gtf2bed_env
from .gtffilter import env as gtffilter_env
from .gtfgenefilter import env as gtfgenefilter_env
from .gunzip import env as gunzip_env
from .hisat2 import env as hisat2_env
from .hmmer import env as hmmer_env
from .htslib import env as htslib_env
from .ilastik import env as ilastik_env
from .interproscan import env as interproscan_env
from .intervals import env as intervals_env
from .kallisto import env as kallisto_env
from .kraken2 import env as kraken2_env
from .krona import env as krona_env
from .macs2 import env as macs2_env
from .mafft import env as mafft_env
from .md5sum import env as md5sum_env
from .megahit import env as megahit_env
from .mindagap import env as mindagap_env
from .minimap2 import env as minimap2_env
from .mmseqs import env as mmseqs_env
from .mosdepth import env as mosdepth_env
from .msisensorpro import env as msisensorpro_env
from .multiqc import env as multiqc_env
from .multiqccustombiotype import env as multiqccustombiotype_env
from .nacho import env as nacho_env
from .nanolyse import env as nanolyse_env
from .nanoplot import env as nanoplot_env
from .ncbigenomedownload import env as ncbigenomedownload_env
from .optitype import env as optitype_env
from .pcgr import env as pcgr_env
from .picard import env as picard_env
from .pigz import env as pigz_env
from .porechop import env as porechop_env
from .preseq import env as preseq_env
from .prodigal import env as prodigal_env
from .prokka import env as prokka_env
from .pydamage import env as pydamage_env
from .qcatch import env as qcatch_env
from .qualimap import env as qualimap_env
from .quarto import env as quarto_env
from .quast import env as quast_env
from .rclone import env as rclone_env
from .rsem import env as rsem_env
from .rseqc import env as rseqc_env
from .rtgtools import env as rtgtools_env
from .s4pred import env as s4pred_env
from .salmon import env as salmon_env
from .sam2lca import env as sam2lca_env
from .samtools import env as samtools_env
from .scanpy import env as scanpy_env
from .seqfu import env as seqfu_env
from .seqkit import env as seqkit_env
from .seqtk import env as seqtk_env
from .simpleaf import env as simpleaf_env
from .snpeff import env as snpeff_env
from .sortmerna import env as sortmerna_env
from .sourcepredict import env as sourcepredict_env
from .sourmash import env as sourmash_env
from .spaceranger import env as spaceranger_env
from .spades import env as spades_env
from .sratools import env as sratools_env
from .star import env as star_env
from .stardist import env as stardist_env
from .stringtie import env as stringtie_env
from .subread import env as subread_env
from .summarizedexperiment import env as summarizedexperiment_env
from .thermorawfileparser import env as thermorawfileparser_env
from .toulligqc import env as toulligqc_env
from .trimgalore import env as trimgalore_env
from .trimmomatic import env as trimmomatic_env
from .trinity import env as trinity_env
from .trust4 import env as trust4_env
from .tx2gene import env as tx2gene_env
from .tximport import env as tximport_env
from .ucsc import env as ucsc_env
from .umicollapse import env as umicollapse_env
from .umitools import env as umitools_env
from .untar import env as untar_env
from .unzip import env as unzip_env
from .vcftools import env as vcftools_env
from .vsearch import env as vsearch_env
from .xz import env as xz_env
from .yara import env as yara_env

env = flyte.TaskEnvironment(
    name="flyte_bio_modules",
    depends_on=[
        amulety_env,
        anndata_env,
        anndatar_env,
        aria2_env,
        art_env,
        bakta_env,
        bbmap_env,
        bcftools_env,
        bedops_env,
        bedtools_env,
        blast_env,
        bowtie2_env,
        bowtie_env,
        bracken_env,
        busco_env,
        bwa_env,
        bwamem2_env,
        cat_env,
        catadditionalfasta_env,
        cellbender_env,
        cellpose_env,
        chromsizes_env,
        clustalo_env,
        cooler_env,
        csvtk_env,
        cutadapt_env,
        damageprofiler_env,
        deeptools_env,
        deepvariant_env,
        deseq2_qc_env,
        diamond_env,
        dragmap_env,
        duckdb_env,
        dupradar_env,
        ensemblvep_env,
        evigene_env,
        fastp_env,
        fastqc_env,
        fastqdl_env,
        fastqutils_env,
        filtlong_env,
        flye_env,
        foldseek_env,
        fq_env,
        gatk4_env,
        gawk_env,
        genomecov_env,
        gffread_env,
        gtf2bed_env,
        gtffilter_env,
        gtfgenefilter_env,
        gunzip_env,
        hisat2_env,
        hmmer_env,
        htslib_env,
        ilastik_env,
        interproscan_env,
        intervals_env,
        kallisto_env,
        kraken2_env,
        krona_env,
        macs2_env,
        mafft_env,
        md5sum_env,
        megahit_env,
        mindagap_env,
        minimap2_env,
        mmseqs_env,
        mosdepth_env,
        msisensorpro_env,
        multiqc_env,
        multiqccustombiotype_env,
        nacho_env,
        nanolyse_env,
        nanoplot_env,
        ncbigenomedownload_env,
        optitype_env,
        pcgr_env,
        picard_env,
        pigz_env,
        porechop_env,
        preseq_env,
        prodigal_env,
        prokka_env,
        pydamage_env,
        qcatch_env,
        qualimap_env,
        quarto_env,
        quast_env,
        rclone_env,
        rsem_env,
        rseqc_env,
        rtgtools_env,
        s4pred_env,
        salmon_env,
        sam2lca_env,
        samtools_env,
        scanpy_env,
        seqfu_env,
        seqkit_env,
        seqtk_env,
        simpleaf_env,
        snpeff_env,
        sortmerna_env,
        sourcepredict_env,
        sourmash_env,
        spades_env,
        sratools_env,
        star_env,
        stardist_env,
        stringtie_env,
        subread_env,
        summarizedexperiment_env,
        thermorawfileparser_env,
        toulligqc_env,
        trimgalore_env,
        trimmomatic_env,
        trinity_env,
        trust4_env,
        tx2gene_env,
        tximport_env,
        ucsc_env,
        umicollapse_env,
        umitools_env,
        untar_env,
        unzip_env,
        vcftools_env,
        vsearch_env,
        xz_env,
        yara_env,
    ],
)

# Tools whose licences restrict use: Cell Ranger and Space Ranger (10x Genomics'
# end-user licence: internal use with 10x products or data, no redistribution)
# and DeepCell's Mesmer (non-commercial use only). They are left out of the
# default `env`; depend on `licensed_env` explicitly once you've checked the terms.
licensed_env = flyte.TaskEnvironment(
    name="flyte_bio_licensed_modules",
    depends_on=[cellranger_env, deepcell_env, spaceranger_env],
)

__all__ = [
    "amulety_env",
    "anndata_env",
    "anndatar_env",
    "aria2_env",
    "art_env",
    "bakta_env",
    "bbmap_env",
    "bcftools_env",
    "bedops_env",
    "bedtools_env",
    "blast_env",
    "bowtie2_env",
    "bowtie_env",
    "bracken_env",
    "busco_env",
    "bwa_env",
    "bwamem2_env",
    "cat_env",
    "catadditionalfasta_env",
    "cellbender_env",
    "cellpose_env",
    "cellranger_env",
    "chromsizes_env",
    "clustalo_env",
    "cooler_env",
    "csvtk_env",
    "cutadapt_env",
    "damageprofiler_env",
    "deepcell_env",
    "deeptools_env",
    "deepvariant_env",
    "deseq2_qc_env",
    "diamond_env",
    "dragmap_env",
    "duckdb_env",
    "dupradar_env",
    "ensemblvep_env",
    "env",
    "evigene_env",
    "fastp_env",
    "fastqc_env",
    "fastqdl_env",
    "fastqutils_env",
    "filtlong_env",
    "flye_env",
    "foldseek_env",
    "fq_env",
    "gatk4_env",
    "gawk_env",
    "genomecov_env",
    "gffread_env",
    "gtf2bed_env",
    "gtffilter_env",
    "gtfgenefilter_env",
    "gunzip_env",
    "hisat2_env",
    "hmmer_env",
    "htslib_env",
    "ilastik_env",
    "interproscan_env",
    "intervals_env",
    "kallisto_env",
    "kraken2_env",
    "krona_env",
    "licensed_env",
    "macs2_env",
    "mafft_env",
    "md5sum_env",
    "megahit_env",
    "mindagap_env",
    "minimap2_env",
    "mmseqs_env",
    "mosdepth_env",
    "msisensorpro_env",
    "multiqc_env",
    "multiqccustombiotype_env",
    "nacho_env",
    "nanolyse_env",
    "nanoplot_env",
    "ncbigenomedownload_env",
    "optitype_env",
    "pcgr_env",
    "picard_env",
    "pigz_env",
    "porechop_env",
    "preseq_env",
    "prodigal_env",
    "prokka_env",
    "pydamage_env",
    "qcatch_env",
    "qualimap_env",
    "quarto_env",
    "quast_env",
    "rclone_env",
    "rsem_env",
    "rseqc_env",
    "rtgtools_env",
    "s4pred_env",
    "salmon_env",
    "sam2lca_env",
    "samtools_env",
    "scanpy_env",
    "seqfu_env",
    "seqkit_env",
    "seqtk_env",
    "simpleaf_env",
    "snpeff_env",
    "sortmerna_env",
    "sourcepredict_env",
    "sourmash_env",
    "spaceranger_env",
    "spades_env",
    "sratools_env",
    "star_env",
    "stardist_env",
    "stringtie_env",
    "subread_env",
    "summarizedexperiment_env",
    "thermorawfileparser_env",
    "toulligqc_env",
    "trimgalore_env",
    "trimmomatic_env",
    "trinity_env",
    "trust4_env",
    "tx2gene_env",
    "tximport_env",
    "ucsc_env",
    "umicollapse_env",
    "umitools_env",
    "untar_env",
    "unzip_env",
    "vcftools_env",
    "vsearch_env",
    "xz_env",
    "yara_env",
]
