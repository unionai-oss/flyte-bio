"""flyte_bio.modules — typed shell-task wrappers around bio CLI tools.

Each submodule wraps one tool family, sharing a single biocontainer image
across its tasks:

- :mod:`flyte_bio.modules.amulety` — AIRR translation and sequence embeddings.
- :mod:`flyte_bio.modules.anndata` — subset an h5ad to listed barcodes.
- :mod:`flyte_bio.modules.anndatar` — h5ad to Seurat / SingleCellExperiment RDS.
- :mod:`flyte_bio.modules.aria2` — file downloads.
- :mod:`flyte_bio.modules.art` — Illumina read simulation.
- :mod:`flyte_bio.modules.assemblyscan` — assembly summary statistics.
- :mod:`flyte_bio.modules.ataqv` — ATAC-seq quality metrics and interactive reports.
- :mod:`flyte_bio.modules.bakta` — bacterial genome annotation (database download, annotate).
- :mod:`flyte_bio.modules.bamtools` — BAM statistics and splitting.
- :mod:`flyte_bio.modules.bbmap` — BBSplit (index build + read binning).
- :mod:`flyte_bio.modules.bcftools` — VCF/BCF utilities (stats, sort, norm, view, filter, annotate, merge, mpileup).
- :mod:`flyte_bio.modules.bedops` — GTF to BED conversion.
- :mod:`flyte_bio.modules.bedtools` — genome arithmetic (intersect, sort,
  merge, genomecov).
- :mod:`flyte_bio.modules.bioawk` — awk with built-in parsing of sequence and interval formats.
- :mod:`flyte_bio.modules.blast` — BLAST databases and blastn search.
- :mod:`flyte_bio.modules.blat` — BLAT alignment (licensed; in ``licensed_env``).
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
- :mod:`flyte_bio.modules.checkqc` — Illumina run-folder QC thresholds.
- :mod:`flyte_bio.modules.chelae` — adapter trimming of FASTQ reads.
- :mod:`flyte_bio.modules.chopper` — long-read quality and length filtering.
- :mod:`flyte_bio.modules.chromap` — fast alignment and preprocessing of chromatin profiles.
- :mod:`flyte_bio.modules.chromsizes` — chromosome sizes and FASTA index.
- :mod:`flyte_bio.modules.clair3` — deep-learning germline small-variant calling from long reads.
- :mod:`flyte_bio.modules.clipkit` — multiple-alignment trimming.
- :mod:`flyte_bio.modules.clustalo` — multiple sequence alignment and guide trees.
- :mod:`flyte_bio.modules.cmaple` — maximum-likelihood phylogenies of closely related sequences.
- :mod:`flyte_bio.modules.collectfeaturecounts` — per-sample featureCounts tables merged into one.
- :mod:`flyte_bio.modules.collectstats` — per-sample read statistics joined into one table.
- :mod:`flyte_bio.modules.cooler` — Hi-C contact matrices.
- :mod:`flyte_bio.modules.crisprcleanr` — correction of gene-independent CRISPR screen effects.
- :mod:`flyte_bio.modules.csvtk` — table joins and concatenation.
- :mod:`flyte_bio.modules.cutadapt` — adapter, primer and quality trimming.
- :mod:`flyte_bio.modules.damageprofiler` — ancient-DNA damage patterns.
- :mod:`flyte_bio.modules.deeptools` — coverage matrices, heatmaps and fingerprint plots.
- :mod:`flyte_bio.modules.deepvariant` — germline small-variant calling.
- :mod:`flyte_bio.modules.dexseq` — differential exon and transcript usage.
- :mod:`flyte_bio.modules.diamond` — protein sequence clustering.
- :mod:`flyte_bio.modules.dragmap` — DRAGMAP short-read alignment (hash table, align).
- :mod:`flyte_bio.modules.duckdb` — CSV/TSV to Parquet conversion.
- :mod:`flyte_bio.modules.eggnogmapper` — functional annotation of proteins with eggNOG.
- :mod:`flyte_bio.modules.ensemblvep` — variant annotation with Ensembl VEP (cache download, annotate).
- :mod:`flyte_bio.modules.epang` — phylogenetic placement on a reference tree.
- :mod:`flyte_bio.modules.evigene` — transcript set reduction (tr2aacds).
- :mod:`flyte_bio.modules.famsa` — fast multiple sequence alignment and guide trees.
- :mod:`flyte_bio.modules.fastqdl` — ENA/SRA FASTQ downloads.
- :mod:`flyte_bio.modules.fastqe` — FASTQ quality summaries as emoji.
- :mod:`flyte_bio.modules.fastqscreen` — screening reads against several genomes.
- :mod:`flyte_bio.modules.fastqutils` — FASTQ validation.
- :mod:`flyte_bio.modules.fibertoolsrs` — m6A and nucleosome calling on Fiber-seq reads.
- :mod:`flyte_bio.modules.filtlong` — long-read quality filtering.
- :mod:`flyte_bio.modules.flye` — long-read de novo assembly.
- :mod:`flyte_bio.modules.foldseek` — protein structure database building and search.
- :mod:`flyte_bio.modules.fq` — FASTQ subsampling and linting.
- :mod:`flyte_bio.modules.gappa` — analysis of phylogenetic placements.
- :mod:`flyte_bio.modules.gatk4` — GATK tools (preprocessing, germline and somatic calling, filtering).
- :mod:`flyte_bio.modules.gawk` — run an awk program over files.
- :mod:`flyte_bio.modules.genmap` — FM index building for mappability.
- :mod:`flyte_bio.modules.gfaffix` — collapsing of shared affixes in variation graphs.
- :mod:`flyte_bio.modules.gffread` — GFF/GTF conversion + transcript FASTA.
- :mod:`flyte_bio.modules.gnusort` — GNU sort of text files.
- :mod:`flyte_bio.modules.graphtyper` — population-scale small-variant genotyping.
- :mod:`flyte_bio.modules.gstama` — TAMA transcript collapse, merge and poly(A) cleanup.
- :mod:`flyte_bio.modules.gtf2bed` — derive a BED12 gene model from a GTF.
- :mod:`flyte_bio.modules.gtffilter` — restrict a GTF to a genome's sequences.
- :mod:`flyte_bio.modules.gtfgenefilter` — keep GTF lines on a genome's sequences (scRNA-seq).
- :mod:`flyte_bio.modules.gunzip` — single-file gzip decompression.
- :mod:`flyte_bio.modules.gzrt` — recover data from corrupt gzipped FASTQs.
- :mod:`flyte_bio.modules.hhsuite` — multiple-alignment format conversion.
- :mod:`flyte_bio.modules.hisat2` — splice-aware short-read alignment (splice sites, build, align).
- :mod:`flyte_bio.modules.hmmer` — profile HMM search (hmmsearch).
- :mod:`flyte_bio.modules.homer` — tag directories, peak calling and peak annotation.
- :mod:`flyte_bio.modules.htslib` — bgzip + tabix.
- :mod:`flyte_bio.modules.ilastik` — pixel classification and multicut segmentation.
- :mod:`flyte_bio.modules.interproscan` — protein domain and family annotation.
- :mod:`flyte_bio.modules.intervals` — scatter-gather genomic intervals.
- :mod:`flyte_bio.modules.iqtree` — maximum-likelihood phylogenies.
- :mod:`flyte_bio.modules.isoseq` — Iso-Seq refine to FLNC reads.
- :mod:`flyte_bio.modules.kallisto` — kallisto index build.
- :mod:`flyte_bio.modules.khmer` — k-mer counting, including unique-k-mer genome size.
- :mod:`flyte_bio.modules.kofamscan` — KEGG Orthology assignment.
- :mod:`flyte_bio.modules.kraken2` — taxonomic classification of reads.
- :mod:`flyte_bio.modules.krona` — interactive taxonomy charts.
- :mod:`flyte_bio.modules.last` — LAST genome alignment, training, splitting, conversion and dot plots.
- :mod:`flyte_bio.modules.leafcutter` — intron clustering for differential splicing.
- :mod:`flyte_bio.modules.lima` — PacBio primer and barcode removal.
- :mod:`flyte_bio.modules.macs2` — ChIP/ATAC-seq peak calling.
- :mod:`flyte_bio.modules.macs3` — model-based peak calling for ChIP-seq and ATAC-seq.
- :mod:`flyte_bio.modules.mafft` — multiple sequence alignment.
- :mod:`flyte_bio.modules.mageck` — CRISPR screen sgRNA counting and gene ranking.
- :mod:`flyte_bio.modules.mat2json` — MATLAB .mat to JSON/CSV conversion.
- :mod:`flyte_bio.modules.md5sum` — MD5 checksums of files.
- :mod:`flyte_bio.modules.medaka` — consensus polishing of nanopore assemblies.
- :mod:`flyte_bio.modules.megahit` — short-read and metagenome assembly.
- :mod:`flyte_bio.modules.metaeuk` — eukaryotic gene prediction in contigs.
- :mod:`flyte_bio.modules.mgnifam` — protein family generation from clustering chunks.
- :mod:`flyte_bio.modules.mindagap` — tile-gap filling and duplicate spot finding.
- :mod:`flyte_bio.modules.minimap2` — long- and short-read alignment (index, align).
- :mod:`flyte_bio.modules.mirdeep2` — read mapping and miRNA discovery.
- :mod:`flyte_bio.modules.mirtop` — miRNA annotation in mirGFF3, with counts and export.
- :mod:`flyte_bio.modules.mirtrace` — quality control for small-RNA sequencing.
- :mod:`flyte_bio.modules.misopy` — MISO annotation indexing.
- :mod:`flyte_bio.modules.mmseqs` — MMseqs2 databases, clustering, taxonomy and tables.
- :mod:`flyte_bio.modules.modkit` — ONT modified-base calling, pileup and repair (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.mosdepth` — read-depth summaries.
- :mod:`flyte_bio.modules.msisensorpro` — microsatellite site scanning.
- :mod:`flyte_bio.modules.multiqc` — aggregate QC report.
- :mod:`flyte_bio.modules.multiqccustombiotype` — biotype counts for MultiQC.
- :mod:`flyte_bio.modules.multiqcsav` — MultiQC reports with Illumina SAV metrics.
- :mod:`flyte_bio.modules.nacho` — NanoString QC and normalisation.
- :mod:`flyte_bio.modules.nanocomp` — long-read dataset comparison.
- :mod:`flyte_bio.modules.nanolyse` — lambda phage read removal.
- :mod:`flyte_bio.modules.nanoplot` — long-read QC plots and statistics.
- :mod:`flyte_bio.modules.ncbigenomedownload` — NCBI genome downloads.
- :mod:`flyte_bio.modules.numorph` — light-sheet preprocessing, resampling and nuclei counting.
- :mod:`flyte_bio.modules.odgi` — building, sorting, layout, drawing and stats of variation graphs.
- :mod:`flyte_bio.modules.optitype` — HLA class I genotyping.
- :mod:`flyte_bio.modules.pairix` — indexing of contact pairs files.
- :mod:`flyte_bio.modules.pairtools` — Hi-C contact pair parsing, sorting, deduplication and filtering.
- :mod:`flyte_bio.modules.pbccs` — PacBio HiFi consensus reads.
- :mod:`flyte_bio.modules.pbjasmine` — 5mC calling on PacBio HiFi reads.
- :mod:`flyte_bio.modules.pbmm2` — minimap2 alignment for PacBio reads.
- :mod:`flyte_bio.modules.pcgr` — PCGR reference bundle download.
- :mod:`flyte_bio.modules.pear` — paired-end read merging (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.phantompeakqualtools` — ChIP-seq strand cross-correlation QC.
- :mod:`flyte_bio.modules.picard` — Picard MarkDuplicates.
- :mod:`flyte_bio.modules.pigz` — parallel gzip compression and decompression.
- :mod:`flyte_bio.modules.pints` — peak calling for nascent transcription data.
- :mod:`flyte_bio.modules.pirate` — prokaryote pangenome gene families.
- :mod:`flyte_bio.modules.porechop` — Nanopore adapter trimming.
- :mod:`flyte_bio.modules.preseq` — library complexity extrapolation.
- :mod:`flyte_bio.modules.prodigal` — prokaryotic gene prediction.
- :mod:`flyte_bio.modules.prokka` — prokaryotic and viral genome annotation.
- :mod:`flyte_bio.modules.pydamage` — per-contig ancient-DNA damage estimates.
- :mod:`flyte_bio.modules.qcatch` — cell-level QC of a simpleaf quantification.
- :mod:`flyte_bio.modules.qualimap` — RNA-seq alignment QC.
- :mod:`flyte_bio.modules.quarto` — render parameterised notebooks.
- :mod:`flyte_bio.modules.quast` — assembly quality assessment.
- :mod:`flyte_bio.modules.racon` — consensus polishing of long-read assemblies.
- :mod:`flyte_bio.modules.rclone` — copy and verify data between storage remotes.
- :mod:`flyte_bio.modules.regtools` — splice junction extraction.
- :mod:`flyte_bio.modules.resolvetaxonomy` — taxonomy tables from taxonomy files or headers.
- :mod:`flyte_bio.modules.ribodetector` — rRNA read removal.
- :mod:`flyte_bio.modules.riker` — one-pass alignment QC metrics.
- :mod:`flyte_bio.modules.rmats` — rMATS prep (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.rsem` — RSEM reference preparation.
- :mod:`flyte_bio.modules.rsemmergecounts` — RSEM results merged into count and TPM matrices.
- :mod:`flyte_bio.modules.rseqc` — RSeQC BAM QC (7 scripts).
- :mod:`flyte_bio.modules.rtgtools` — RTG Tools SDF conversion.
- :mod:`flyte_bio.modules.rundbcan` — CAZyme annotation and its database.
- :mod:`flyte_bio.modules.rustqc` — one-pass RNA-seq alignment QC.
- :mod:`flyte_bio.modules.s4pred` — protein secondary structure prediction.
- :mod:`flyte_bio.modules.salmon` — transcript quantification (index, quant).
- :mod:`flyte_bio.modules.sam2lca` — lowest-common-ancestor taxonomic assignment.
- :mod:`flyte_bio.modules.scanpy` — AnnData conversions (simpleaf to h5ad, concatenation).
- :mod:`flyte_bio.modules.simpleaf` — single-cell quantification with alevin-fry (index, quant).
- :mod:`flyte_bio.modules.samtools` — alignment sort/index/stats/faidx.
- :mod:`flyte_bio.modules.seacr` — CUT&RUN/CUT&Tag peak calling.
- :mod:`flyte_bio.modules.seq2hla` — HLA typing and expression from RNA-seq.
- :mod:`flyte_bio.modules.seqcluster` — collapsing of small-RNA reads.
- :mod:`flyte_bio.modules.seqfu` — sequence-length statistics.
- :mod:`flyte_bio.modules.seqkit` — FASTA/FASTQ statistics and pattern replacement.
- :mod:`flyte_bio.modules.seqtk` — FASTQ subsampling.
- :mod:`flyte_bio.modules.sequali` — sequencing-read QC reports.
- :mod:`flyte_bio.modules.seqwish` — variation graph induction from alignments.
- :mod:`flyte_bio.modules.smoothxg` — variation graph smoothing.
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
- :mod:`flyte_bio.modules.suppa` — alternative splicing events and differential splicing.
- :mod:`flyte_bio.modules.sylph` — metagenomic profiling by k-mer containment ANI.
- :mod:`flyte_bio.modules.sylphtax` — taxonomic profiles from sylph results.
- :mod:`flyte_bio.modules.taxonkit` — NCBI taxonomy lineages.
- :mod:`flyte_bio.modules.thermorawfileparser` — Thermo RAW conversion.
- :mod:`flyte_bio.modules.toulligqc` — Nanopore run QC reports.
- :mod:`flyte_bio.modules.transdecoder` — coding regions in transcripts.
- :mod:`flyte_bio.modules.trimgalore` — adapter/quality trimming.
- :mod:`flyte_bio.modules.trimmomatic` — Illumina adapter and quality trimming.
- :mod:`flyte_bio.modules.trinity` — de novo transcriptome assembly.
- :mod:`flyte_bio.modules.trust4` — TCR/BCR repertoire reconstruction.
- :mod:`flyte_bio.modules.tx2gene` — transcript → gene table from a GTF.
- :mod:`flyte_bio.modules.tximport` — count/TPM matrices from quantifications.
- :mod:`flyte_bio.modules.ucsc` — bedClip / bedGraphToBigWig.
- :mod:`flyte_bio.modules.ultra` — annotation-guided long-read splice alignment.
- :mod:`flyte_bio.modules.umicollapse` — UMI deduplication.
- :mod:`flyte_bio.modules.umitools` — UMI extraction and deduplication.
- :mod:`flyte_bio.modules.untar` — tar archive extraction.
- :mod:`flyte_bio.modules.untarfiles` — tar extraction with a file listing.
- :mod:`flyte_bio.modules.unzip` — zip archive extraction.
- :mod:`flyte_bio.modules.unzipfiles` — archive extraction with a file list.
- :mod:`flyte_bio.modules.vcflib` — VCF filtering by field expressions.
- :mod:`flyte_bio.modules.vcftools` — VCF summaries (Ts/Tv, filters).
- :mod:`flyte_bio.modules.vsearch` — sequence clustering.
- :mod:`flyte_bio.modules.wfmash` — all-to-all genome alignment for pangenome graphs.
- :mod:`flyte_bio.modules.wget` — file download with GNU Wget.
- :mod:`flyte_bio.modules.wipertools` — repair malformed FASTQs in parallel chunks.
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
from .assemblyscan import env as assemblyscan_env
from .ataqv import env as ataqv_env
from .bakta import env as bakta_env
from .bamtools import env as bamtools_env
from .bbmap import env as bbmap_env
from .bcftools import env as bcftools_env
from .bedops import env as bedops_env
from .bedtools import env as bedtools_env
from .bedtools import genomecov_env
from .bioawk import env as bioawk_env
from .blast import env as blast_env
from .blat import env as blat_env
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
from .checkqc import env as checkqc_env
from .chelae import env as chelae_env
from .chopper import env as chopper_env
from .chromap import env as chromap_env
from .chromsizes import env as chromsizes_env
from .clair3 import env as clair3_env
from .clipkit import env as clipkit_env
from .clustalo import env as clustalo_env
from .cmaple import env as cmaple_env
from .collectfeaturecounts import env as collectfeaturecounts_env
from .collectstats import env as collectstats_env
from .cooler import env as cooler_env
from .crisprcleanr import env as crisprcleanr_env
from .csvtk import env as csvtk_env
from .cutadapt import env as cutadapt_env
from .damageprofiler import env as damageprofiler_env
from .deepcell import env as deepcell_env
from .deeptools import env as deeptools_env
from .deepvariant import env as deepvariant_env
from .deseq2_qc import env as deseq2_qc_env
from .dexseq import env as dexseq_env
from .diamond import env as diamond_env
from .dragmap import env as dragmap_env
from .duckdb import env as duckdb_env
from .dupradar import env as dupradar_env
from .eggnogmapper import env as eggnogmapper_env
from .ensemblvep import env as ensemblvep_env
from .epang import env as epang_env
from .evigene import env as evigene_env
from .famsa import env as famsa_env
from .fastp import env as fastp_env
from .fastqc import env as fastqc_env
from .fastqdl import env as fastqdl_env
from .fastqe import env as fastqe_env
from .fastqscreen import env as fastqscreen_env
from .fastqutils import env as fastqutils_env
from .fibertoolsrs import env as fibertoolsrs_env
from .filtlong import env as filtlong_env
from .flye import env as flye_env
from .foldseek import env as foldseek_env
from .fq import env as fq_env
from .gappa import env as gappa_env
from .gatk4 import env as gatk4_env
from .gawk import env as gawk_env
from .genmap import env as genmap_env
from .gfaffix import env as gfaffix_env
from .gffread import env as gffread_env
from .gnusort import env as gnusort_env
from .graphtyper import env as graphtyper_env
from .gstama import env as gstama_env
from .gtf2bed import env as gtf2bed_env
from .gtffilter import env as gtffilter_env
from .gtfgenefilter import env as gtfgenefilter_env
from .gunzip import env as gunzip_env
from .gzrt import env as gzrt_env
from .hhsuite import env as hhsuite_env
from .hisat2 import env as hisat2_env
from .hmmer import env as hmmer_env
from .homer import env as homer_env
from .htslib import env as htslib_env
from .ilastik import env as ilastik_env
from .interproscan import env as interproscan_env
from .intervals import env as intervals_env
from .iqtree import env as iqtree_env
from .isoseq import env as isoseq_env
from .kallisto import env as kallisto_env
from .khmer import env as khmer_env
from .kofamscan import env as kofamscan_env
from .kraken2 import env as kraken2_env
from .krona import env as krona_env
from .last import env as last_env
from .leafcutter import env as leafcutter_env
from .lima import env as lima_env
from .macs2 import env as macs2_env
from .macs3 import env as macs3_env
from .mafft import env as mafft_env
from .mageck import env as mageck_env
from .mat2json import env as mat2json_env
from .md5sum import env as md5sum_env
from .medaka import env as medaka_env
from .megahit import env as megahit_env
from .metaeuk import env as metaeuk_env
from .mgnifam import env as mgnifam_env
from .mindagap import env as mindagap_env
from .minimap2 import env as minimap2_env
from .mirdeep2 import env as mirdeep2_env
from .mirtop import env as mirtop_env
from .mirtrace import env as mirtrace_env
from .misopy import env as misopy_env
from .mmseqs import env as mmseqs_env
from .modkit import env as modkit_env
from .mosdepth import env as mosdepth_env
from .msisensorpro import env as msisensorpro_env
from .multiqc import env as multiqc_env
from .multiqccustombiotype import env as multiqccustombiotype_env
from .multiqcsav import env as multiqcsav_env
from .nacho import env as nacho_env
from .nanocomp import env as nanocomp_env
from .nanolyse import env as nanolyse_env
from .nanoplot import env as nanoplot_env
from .ncbigenomedownload import env as ncbigenomedownload_env
from .numorph import env as numorph_env
from .odgi import env as odgi_env
from .optitype import env as optitype_env
from .pairix import env as pairix_env
from .pairtools import env as pairtools_env
from .pbccs import env as pbccs_env
from .pbjasmine import env as pbjasmine_env
from .pbmm2 import env as pbmm2_env
from .pcgr import env as pcgr_env
from .pear import env as pear_env
from .phantompeakqualtools import env as phantompeakqualtools_env
from .picard import env as picard_env
from .pigz import env as pigz_env
from .pints import env as pints_env
from .pirate import env as pirate_env
from .porechop import env as porechop_env
from .preseq import env as preseq_env
from .prodigal import env as prodigal_env
from .prokka import env as prokka_env
from .pydamage import env as pydamage_env
from .qcatch import env as qcatch_env
from .qualimap import env as qualimap_env
from .quarto import env as quarto_env
from .quast import env as quast_env
from .racon import env as racon_env
from .rclone import env as rclone_env
from .regtools import env as regtools_env
from .resolvetaxonomy import env as resolvetaxonomy_env
from .ribodetector import env as ribodetector_env
from .riker import env as riker_env
from .rmats import env as rmats_env
from .rsem import env as rsem_env
from .rsemmergecounts import env as rsemmergecounts_env
from .rseqc import env as rseqc_env
from .rtgtools import env as rtgtools_env
from .rundbcan import env as rundbcan_env
from .rustqc import env as rustqc_env
from .s4pred import env as s4pred_env
from .salmon import env as salmon_env
from .sam2lca import env as sam2lca_env
from .samtools import env as samtools_env
from .scanpy import env as scanpy_env
from .seacr import env as seacr_env
from .seq2hla import env as seq2hla_env
from .seqcluster import env as seqcluster_env
from .seqfu import env as seqfu_env
from .seqkit import env as seqkit_env
from .seqtk import env as seqtk_env
from .sequali import env as sequali_env
from .seqwish import env as seqwish_env
from .simpleaf import env as simpleaf_env
from .smoothxg import env as smoothxg_env
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
from .suppa import env as suppa_env
from .sylph import env as sylph_env
from .sylphtax import env as sylphtax_env
from .taxonkit import env as taxonkit_env
from .thermorawfileparser import env as thermorawfileparser_env
from .toulligqc import env as toulligqc_env
from .transdecoder import env as transdecoder_env
from .trimgalore import env as trimgalore_env
from .trimmomatic import env as trimmomatic_env
from .trinity import env as trinity_env
from .trust4 import env as trust4_env
from .tx2gene import env as tx2gene_env
from .tximport import env as tximport_env
from .ucsc import env as ucsc_env
from .ultra import env as ultra_env
from .umicollapse import env as umicollapse_env
from .umitools import env as umitools_env
from .untar import env as untar_env
from .untarfiles import env as untarfiles_env
from .unzip import env as unzip_env
from .unzipfiles import env as unzipfiles_env
from .vcflib import env as vcflib_env
from .vcftools import env as vcftools_env
from .vsearch import env as vsearch_env
from .wfmash import env as wfmash_env
from .wget import env as wget_env
from .wipertools import env as wipertools_env
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
        assemblyscan_env,
        ataqv_env,
        bakta_env,
        bamtools_env,
        bbmap_env,
        bcftools_env,
        bedops_env,
        bedtools_env,
        bioawk_env,
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
        checkqc_env,
        chelae_env,
        chopper_env,
        chromap_env,
        chromsizes_env,
        clair3_env,
        clipkit_env,
        clustalo_env,
        cmaple_env,
        collectfeaturecounts_env,
        collectstats_env,
        cooler_env,
        crisprcleanr_env,
        csvtk_env,
        cutadapt_env,
        damageprofiler_env,
        deeptools_env,
        deepvariant_env,
        deseq2_qc_env,
        dexseq_env,
        diamond_env,
        dragmap_env,
        duckdb_env,
        dupradar_env,
        eggnogmapper_env,
        ensemblvep_env,
        epang_env,
        evigene_env,
        famsa_env,
        fastp_env,
        fastqc_env,
        fastqdl_env,
        fastqe_env,
        fastqscreen_env,
        fastqutils_env,
        fibertoolsrs_env,
        filtlong_env,
        flye_env,
        foldseek_env,
        fq_env,
        gappa_env,
        gatk4_env,
        gawk_env,
        genmap_env,
        genomecov_env,
        gfaffix_env,
        gffread_env,
        gnusort_env,
        graphtyper_env,
        gstama_env,
        gtf2bed_env,
        gtffilter_env,
        gtfgenefilter_env,
        gunzip_env,
        gzrt_env,
        hhsuite_env,
        hisat2_env,
        hmmer_env,
        homer_env,
        htslib_env,
        ilastik_env,
        interproscan_env,
        intervals_env,
        iqtree_env,
        isoseq_env,
        kallisto_env,
        khmer_env,
        kofamscan_env,
        kraken2_env,
        krona_env,
        last_env,
        leafcutter_env,
        lima_env,
        macs2_env,
        macs3_env,
        mafft_env,
        mageck_env,
        mat2json_env,
        md5sum_env,
        medaka_env,
        megahit_env,
        metaeuk_env,
        mgnifam_env,
        mindagap_env,
        minimap2_env,
        mirdeep2_env,
        mirtop_env,
        mirtrace_env,
        misopy_env,
        mmseqs_env,
        mosdepth_env,
        msisensorpro_env,
        multiqc_env,
        multiqccustombiotype_env,
        multiqcsav_env,
        nacho_env,
        nanocomp_env,
        nanolyse_env,
        nanoplot_env,
        ncbigenomedownload_env,
        numorph_env,
        odgi_env,
        optitype_env,
        pairix_env,
        pairtools_env,
        pbccs_env,
        pbjasmine_env,
        pbmm2_env,
        pcgr_env,
        phantompeakqualtools_env,
        picard_env,
        pigz_env,
        pints_env,
        pirate_env,
        porechop_env,
        preseq_env,
        prodigal_env,
        prokka_env,
        pydamage_env,
        qcatch_env,
        qualimap_env,
        quarto_env,
        quast_env,
        racon_env,
        rclone_env,
        regtools_env,
        resolvetaxonomy_env,
        ribodetector_env,
        riker_env,
        rsem_env,
        rsemmergecounts_env,
        rseqc_env,
        rtgtools_env,
        rundbcan_env,
        rustqc_env,
        s4pred_env,
        salmon_env,
        sam2lca_env,
        samtools_env,
        scanpy_env,
        seacr_env,
        seq2hla_env,
        seqcluster_env,
        seqfu_env,
        seqkit_env,
        seqtk_env,
        sequali_env,
        seqwish_env,
        simpleaf_env,
        smoothxg_env,
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
        suppa_env,
        sylph_env,
        sylphtax_env,
        taxonkit_env,
        thermorawfileparser_env,
        toulligqc_env,
        transdecoder_env,
        trimgalore_env,
        trimmomatic_env,
        trinity_env,
        trust4_env,
        tx2gene_env,
        tximport_env,
        ucsc_env,
        ultra_env,
        umicollapse_env,
        umitools_env,
        untar_env,
        untarfiles_env,
        unzip_env,
        unzipfiles_env,
        vcflib_env,
        vcftools_env,
        vsearch_env,
        wfmash_env,
        wget_env,
        wipertools_env,
        xz_env,
        yara_env,
    ],
)

# Tools whose licences restrict use, left out of the default `env`:
# - Cell Ranger and Space Ranger: 10x Genomics' end-user licence (internal use
#   with 10x products or data, no redistribution).
# - DeepCell's Mesmer, PEAR and rMATS: non-commercial use only.
# - BLAT: free for academic, non-profit and personal use; commercial use needs a
#   Kent Informatics licence.
# - modkit: Oxford Nanopore's public licence (research use).
# Depend on `licensed_env` explicitly once you've checked the terms.
licensed_env = flyte.TaskEnvironment(
    name="flyte_bio_licensed_modules",
    depends_on=[blat_env, cellranger_env, deepcell_env, modkit_env, pear_env, rmats_env, spaceranger_env],
)

__all__ = [
    "amulety_env",
    "anndata_env",
    "anndatar_env",
    "aria2_env",
    "art_env",
    "assemblyscan_env",
    "ataqv_env",
    "bakta_env",
    "bamtools_env",
    "bbmap_env",
    "bcftools_env",
    "bedops_env",
    "bedtools_env",
    "bioawk_env",
    "blast_env",
    "blat_env",
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
    "checkqc_env",
    "chelae_env",
    "chopper_env",
    "chromap_env",
    "chromsizes_env",
    "clair3_env",
    "clipkit_env",
    "clustalo_env",
    "cmaple_env",
    "collectfeaturecounts_env",
    "collectstats_env",
    "cooler_env",
    "crisprcleanr_env",
    "csvtk_env",
    "cutadapt_env",
    "damageprofiler_env",
    "deepcell_env",
    "deeptools_env",
    "deepvariant_env",
    "deseq2_qc_env",
    "dexseq_env",
    "diamond_env",
    "dragmap_env",
    "duckdb_env",
    "dupradar_env",
    "eggnogmapper_env",
    "ensemblvep_env",
    "env",
    "epang_env",
    "evigene_env",
    "famsa_env",
    "fastp_env",
    "fastqc_env",
    "fastqdl_env",
    "fastqe_env",
    "fastqscreen_env",
    "fastqutils_env",
    "fibertoolsrs_env",
    "filtlong_env",
    "flye_env",
    "foldseek_env",
    "fq_env",
    "gappa_env",
    "gatk4_env",
    "gawk_env",
    "genmap_env",
    "genomecov_env",
    "gfaffix_env",
    "gffread_env",
    "gnusort_env",
    "graphtyper_env",
    "gstama_env",
    "gtf2bed_env",
    "gtffilter_env",
    "gtfgenefilter_env",
    "gunzip_env",
    "gzrt_env",
    "hhsuite_env",
    "hisat2_env",
    "hmmer_env",
    "homer_env",
    "htslib_env",
    "ilastik_env",
    "interproscan_env",
    "intervals_env",
    "iqtree_env",
    "isoseq_env",
    "kallisto_env",
    "khmer_env",
    "kofamscan_env",
    "kraken2_env",
    "krona_env",
    "last_env",
    "leafcutter_env",
    "licensed_env",
    "lima_env",
    "macs2_env",
    "macs3_env",
    "mafft_env",
    "mageck_env",
    "mat2json_env",
    "md5sum_env",
    "medaka_env",
    "megahit_env",
    "metaeuk_env",
    "mgnifam_env",
    "mindagap_env",
    "minimap2_env",
    "mirdeep2_env",
    "mirtop_env",
    "mirtrace_env",
    "misopy_env",
    "mmseqs_env",
    "modkit_env",
    "mosdepth_env",
    "msisensorpro_env",
    "multiqc_env",
    "multiqccustombiotype_env",
    "multiqcsav_env",
    "nacho_env",
    "nanocomp_env",
    "nanolyse_env",
    "nanoplot_env",
    "ncbigenomedownload_env",
    "numorph_env",
    "odgi_env",
    "optitype_env",
    "pairix_env",
    "pairtools_env",
    "pbccs_env",
    "pbjasmine_env",
    "pbmm2_env",
    "pcgr_env",
    "pear_env",
    "phantompeakqualtools_env",
    "picard_env",
    "pigz_env",
    "pints_env",
    "pirate_env",
    "porechop_env",
    "preseq_env",
    "prodigal_env",
    "prokka_env",
    "pydamage_env",
    "qcatch_env",
    "qualimap_env",
    "quarto_env",
    "quast_env",
    "racon_env",
    "rclone_env",
    "regtools_env",
    "resolvetaxonomy_env",
    "ribodetector_env",
    "riker_env",
    "rmats_env",
    "rsem_env",
    "rsemmergecounts_env",
    "rseqc_env",
    "rtgtools_env",
    "rundbcan_env",
    "rustqc_env",
    "s4pred_env",
    "salmon_env",
    "sam2lca_env",
    "samtools_env",
    "scanpy_env",
    "seacr_env",
    "seq2hla_env",
    "seqcluster_env",
    "seqfu_env",
    "seqkit_env",
    "seqtk_env",
    "sequali_env",
    "seqwish_env",
    "simpleaf_env",
    "smoothxg_env",
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
    "suppa_env",
    "sylph_env",
    "sylphtax_env",
    "taxonkit_env",
    "thermorawfileparser_env",
    "toulligqc_env",
    "transdecoder_env",
    "trimgalore_env",
    "trimmomatic_env",
    "trinity_env",
    "trust4_env",
    "tx2gene_env",
    "tximport_env",
    "ucsc_env",
    "ultra_env",
    "umicollapse_env",
    "umitools_env",
    "untar_env",
    "untarfiles_env",
    "unzip_env",
    "unzipfiles_env",
    "vcflib_env",
    "vcftools_env",
    "vsearch_env",
    "wfmash_env",
    "wget_env",
    "wipertools_env",
    "xz_env",
    "yara_env",
]
