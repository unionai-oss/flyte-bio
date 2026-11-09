# Modules

`flyte_bio.modules` wraps bioinformatics tools as Flyte tasks, one Python file
per tool family.

## How a module works

Each tool is a Flyte shell task (`flyte.extras.shell`) that runs the tool's
command line in the container upstream pins for it. Usually that's a
[biocontainer](https://biocontainers.pro/) or a Seqera Wave image. Your own
image never needs the tools installed.

A module exposes:

- **Task functions** you `await` from your own task, e.g.
  `await samtools_sort(bam=bam)`. Simple tools are the shell task itself.
  Tools with several outputs get a thin async wrapper that returns a small
  dataclass (e.g. `MarkDuplicatesResult`) or picks files out of the tool's
  output directory by name.
- **An `env`** (a `flyte.TaskEnvironment`) holding the module's tasks. A
  family whose tools need different containers has one environment per
  container and an aggregate `env`. `flyte_bio.modules.env` (re-exported as
  `flyte_bio.env`) depends on every module's `env`, so one `depends_on` entry
  brings in everything.

Conventions shared by every module:

- **Files keep their names.** Input files are staged under their original
  names, and tools that care about names (indexes next to their files, tools
  that record names in their outputs) get the inputs they expect.
- **Outputs keep upstream's names.** When a tool writes several files, they
  land in one output directory under the names the upstream module gives them
  (e.g. `<prefix>.md.cram`, `<prefix>.md.cram.metrics`).
- **`args` passes extra options straight to the tool**, as upstream's
  `ext.args` does. It's split on whitespace, so a single option value can't
  contain spaces.
- **Resources are sized for typical inputs.** Each module sets its tasks'
  CPU and memory near the top of its file, e.g. `bwa.MEM_RESOURCES`. A
  tool's thread count and memory options (`-t`, `-Xmx`, sort buffers) are
  read from the container's limits when it runs (`flyte_bio.limits`), as
  upstream's are from `task.cpus` and `task.memory`. So a tool that runs with
  overridden resources uses them. A tool's task can be overridden
  (`bwa.bwa_mem_cmd.as_task().override(resources=...)`), but the module
  functions don't take resources yet.
- **Helper scripts are vendored unchanged.** Python, R and Perl helpers that a
  container doesn't include live in `flyte_bio.scripts` and are passed to the
  task as files.

Each module's tests check its outputs against the upstream module's own test
cases (md5 snapshots wherever the output is reproducible). See
[contributing](contributing.md).

## Available modules

Pipeline columns show which pipelines use each module.

| Module | What it does | Functions | Used by |
|---|---|---|---|
| `amulety` | Translate AIRR rearrangements with IgBLAST; embed BCR/TCR sequences with protein language models | `amulety_translate`, `amulety_embed` | — |
| `anndata` | Keep an h5ad's cells with listed barcodes | `anndata_barcodes` | scrnaseq |
| `anndatar` | h5ad to Seurat and SingleCellExperiment RDS | `anndatar_convert` | scrnaseq |
| `aria2` | Download a file over HTTP(S) or FTP | `aria2` | — |
| `art` | Simulate Illumina reads from a reference | `art_illumina` | — |
| `ashlar` | Stitch and register multiplexed tissue image cycles into one pyramidal OME-TIFF, with optional illumination correction | `ashlar` | — |
| `assemblyscan` | Assembly summary statistics (contig counts, lengths, N50, GC) as TSV or JSON | `assemblyscan` | — |
| `ataqv` | ATAC-seq quality metrics per sample, and the interactive HTML report across samples | `ataqv_ataqv`, `ataqv_mkarv` | — |
| `autocycler` | Consensus long-read bacterial assembly: subsample reads, compress, cluster, trim, resolve and combine assemblies | `autocycler_subsample`, `autocycler_compress`, `autocycler_cluster`, `autocycler_trim`, `autocycler_resolve`, `autocycler_combine` | — |
| `backsub` | Channel-wise background subtraction of a registered multiplexed image, driven by a markers CSV | `backsub` | — |
| `bakta` | Annotate a bacterial genome, plasmid or MAG; download its database | `bakta_db_download`, `bakta` | — |
| `bamtools` | BAM statistics; split BAMs by reference, mapping, pairing or tag, as a `.tar.gz` of BAMs | `bamtools_stats`, `bamtools_split`, `bamtools_convert` | — |
| `bases2fastq` | Element AVITI run-folder demultiplexing to FASTQ (Element EULA) | `bases2fastq` | — |
| `basicpy` | Estimate dark-field and flat-field illumination profiles of a tiled microscopy image (BaSiC) | `basicpy` | — |
| `bbmap` | BBTools: BBMap indexing and alignment; BBSplit binning of reads by reference; BBDuk trimming; BBNorm read-depth normalisation; Clumpify read grouping; Pileup coverage stats; FilterByName read selection; Repair re-pairing of out-of-order mates | `bbsplit_index`, `bbsplit`, `bbduk`, `bbnorm`, `bbmap_filterbyname`, `bbmap_repair`, `bbmap_index`, `bbmap_align`, `bbmap_pileup`, `bbmap_clumpify`, `bbmap_sendsketch` | rnaseq |
| `bcftools` | VCF/BCF utilities: statistics, sort, normalise, view, filter, annotate, merge, mpileup calling, index, concatenate, reheader, consensus sequence, set genotypes | `bcftools_stats`, `bcftools_sort`, `bcftools_norm`, `bcftools_view`, `bcftools_filter`, `bcftools_annotate`, `bcftools_merge`, `bcftools_mpileup`, `bcftools_index`, `bcftools_concat`, `bcftools_reheader`, `bcftools_consensus`, `bcftools_pluginsetgt` | variant_calling |
| `bcl2fastq` | Illumina BCL run-folder demultiplexing to FASTQ (Illumina EULA) | `bcl2fastq` | — |
| `bclconvert` | Illumina BCL Convert run-folder demultiplexing to FASTQ (Illumina EULA) | `bclconvert` | — |
| `bedops` | GTF to BED conversion | `bedops_gtf2bed` | — |
| `bedtools` | Genome arithmetic; extract or mask a BED's intervals in a FASTA; convert a BAM to BED; complement a BED against a genome | `bedtools_intersect`, `bedtools_sort`, `bedtools_merge`, `bedtools_genomecov`, `bedtools_getfasta`, `bedtools_maskfasta`, `bedtools_bamtobed`, `bedtools_complement` | rnaseq |
| `bff` | Demultiplex cell-hashing (HTO) counts with cellhashR's BFF callers | `bff` | — |
| `bftools` | Bio-Formats tools: dump an image's OME-XML metadata | `bftools_showinf` | — |
| `bioawk` | awk with built-in parsing of FASTA, FASTQ, BED, SAM, VCF and GFF records | `bioawk` | — |
| `blast` | Build BLAST databases and search nucleotide sequences with blastn | `blast_makeblastdb`, `blast_blastn` | — |
| `blat` | BLAT alignment of sequences against a database. Academic, non-profit and personal use only: in `licensed_env`, not the default `env` | `blat` | — |
| `bowtie` | Short-read alignment with Bowtie 1: index build, alignment to BAM | `bowtie1_build`, `bowtie1_align` | — |
| `bowtie2` | Bowtie 2 short-read alignment: build an index, align reads through samtools to BAM/SAM/CRAM, with the log and optional unaligned reads | `bowtie2_build`, `bowtie2_align` | — |
| `bracken` | Abundance re-estimation from a Kraken2 report | `bracken` | — |
| `busco` | Assembly, proteome and transcriptome completeness against a BUSCO lineage; download lineages | `busco_download`, `busco` | — |
| `bwa` | BWA short-read alignment: index, BWA-MEM, and BWA-backtrack (`aln`, then `samse`/`sampe` to sorted BAM) | `bwa_index`, `bwa_mem`, `bwa_aln`, `bwa_samse`, `bwa_sampe` | variant_calling |
| `bwamem2` | BWA-MEM2 index and short-read alignment, piped through samtools sort or view | `bwamem2_index`, `bwamem2_mem` | — |
| `canu` | De novo assembly of long reads (PacBio CLR, HiFi, Nanopore) | `canu` | — |
| `cat` | Concatenate a sample's FASTQs; concatenate any files, (de)compressing by extension | `cat_fastq`, `cat_cat`, `concatenate` | rnaseq |
| `catadditionalfasta` | Append extra sequences (e.g. spike-ins) to a genome and GTF | `cat_additional_fasta` | rnaseq |
| `cellbender` | Remove ambient RNA and empty droplets (CPU or GPU) | `cellbender_removebackground` | scrnaseq |
| `cellpose` | Cell segmentation of microscopy images | `cellpose` | — |
| `cellranger` | 10x Cell Ranger V(D)J reference building and immune profiling. Licensed: in `licensed_env`, not the default `env` | `cellranger_mkvdjref`, `cellranger_vdj` | — |
| `cellsnp` | Per-cell allele counts at known SNPs (cellsnp-lite mode 1a) | `cellsnp_modea` | — |
| `checkqc` | Check an Illumina run folder (archived) against QC thresholds | `checkqc` | — |
| `chelae` | Adapter trimming of FASTQ reads | `chelae_trim` | — |
| `chopper` | Quality, length and contaminant filtering of long reads | `chopper` | — |
| `chromap` | Fast alignment and preprocessing of ChIP-seq, ATAC-seq and Hi-C reads; index building | `chromap_index`, `chromap_chromap` | — |
| `chromsizes` | Chromosome sizes and FASTA index (`.fai`, `.gzi` if bgzipped) | `chrom_sizes` | — |
| `clair3` | Deep-learning germline small-variant calling from long reads | `clair3` | — |
| `clipkit` | Trim multiple sequence alignments (ClipKIT), with per-site log, complement and JSON report | `clipkit` | — |
| `clustalo` | Multiple sequence alignment and guide trees with Clustal Omega | `clustalo_align`, `clustalo_guidetree` | — |
| `cmaple` | Maximum-likelihood phylogeny of closely related sequences (CMAPLE), optionally from a starting tree | `cmaple` | — |
| `cnvkit` | Copy-number calling from targeted and whole-genome sequencing (CNVkit): target BEDs and the batch pipeline (references, panels of normals, copy ratios and segments) | `cnvkit_target`, `cnvkit_batch` | — |
| `collectfeaturecounts` | Collect per-sample featureCounts tables into one long counts table with per-sample TPMs | `collect_featurecounts` | — |
| `collectstats` | Join per-sample read statistics (trimming, BBDuk, idxstats, feature counts, merge table) into one table | `collect_stats` | — |
| `cooler` | Hi-C contact matrices: bin pairs into a .cool, make genomic bins, balance, dump tables, zoomify to a multi-resolution .mcool | `cooler_cload`, `cooler_digest`, `cooler_makebins`, `cooler_balance`, `cooler_dump`, `cooler_zoomify` | — |
| `coreograph` | Detect and crop tissue-microarray cores (U-Net), with masks, TMA map and centroids | `coreograph` | — |
| `crisprcleanr` | Correction of gene-independent copy-number effects in CRISPR screen counts | `crisprcleanr_normalize` | — |
| `csvtk` | Join tables on key fields; concatenate tables, converting the delimiter | `csvtk_join`, `csvtk_concat` | — |
| `cutadapt` | Adapter, primer and quality trimming of FASTQ reads (single- or paired-end) | `cutadapt` | — |
| `damageprofiler` | Ancient-DNA damage patterns (misincorporation, read lengths) of a BAM | `damageprofiler` | — |
| `deepcell` | Whole-cell and nuclear segmentation with Mesmer. Non-commercial licence: in `licensed_env`, not the default `env` | `deepcell_mesmer` | — |
| `deeptools` | bigWig coverage tracks, profile plots, coverage matrices, heatmaps and fingerprint plots for ChIP/ATAC-seq; BAM bin summaries with correlation and PCA plots | `deeptools_computematrix`, `deeptools_plotheatmap`, `deeptools_plotprofile`, `deeptools_plotfingerprint`, `deeptools_multibamsummary`, `deeptools_plotcorrelation`, `deeptools_plotpca`, `deeptools_bamcoverage` | — |
| `deepvariant` | Germline small-variant calling with DeepVariant (VCF and gVCF) | `deepvariant_rundeepvariant` | — |
| `demuxem` | Demultiplex cell-hashing data with demuxEM | `demuxem` | — |
| `deseq2_qc` | Sample-level QC of a count matrix with DESeq2 (PCA, distances) | `deseq2_qc` | rnaseq |
| `dexseq` | Differential exon and transcript usage with DEXSeq: annotation flattening, exon counting, DEU and DTU tests | `dexseq_prepareannotation`, `dexseq_count`, `dexseq_exon`, `dexseq_dtu` | — |
| `diamond` | Cluster protein sequences (FASTA or DIAMOND database) | `diamond_cluster`, `diamond_makedb`, `diamond_blastp`, `diamond_blastx` | — |
| `dragmap` | Short-read alignment with DRAGMAP: hash table build, alignment to BAM/CRAM | `dragmap_hashtable`, `dragmap_align` | — |
| `dragonflye` | Long-read bacterial assembly (Flye, Miniasm or Raven) with long- and short-read polishing | `dragonflye` | — |
| `duckdb` | Convert a CSV/TSV (optionally compressed) to Parquet | `duckdb_table2parquet` | — |
| `dupradar` | Duplication rate vs. expression for RNA-seq BAMs | `dupradar` | rnaseq |
| `eggnogmapper` | Functional annotation of proteins (orthologs, GO, EC, KEGG, CAZy, PFAMs) with eggNOG-mapper | `eggnogmapper` | — |
| `ensemblvep` | Variant annotation with Ensembl VEP, with its cache download | `ensemblvep_download`, `ensemblvep_vep` | — |
| `epang` | Phylogenetic placement of query sequences on a reference tree, and splitting reference/query alignments (EPA-ng) | `epang_place`, `epang_split` | — |
| `evigene` | Reduce an over-assembled transcriptome to non-redundant transcripts, CDS and proteins | `evigene_tr2aacds` | — |
| `falco` | Read QC reports (FastQC drop-in) | `falco` | — |
| `famsa` | Fast multiple sequence alignment (FAMSA), optionally from a guide tree; export FAMSA's guide tree | `famsa_align`, `famsa_guidetree` | — |
| `fastp` | Adapter and quality trimming with QC reports (single-end, paired-end, interleaved) | `fastp` | — |
| `fastqc` | Read quality reports | `fastqc`, `fastqc_reads` | rnaseq, scrnaseq, variant_calling |
| `fastqdl` | Download FASTQ files and run metadata from ENA or SRA | `fastqdl` | — |
| `fastqe` | FASTQ quality summaries as emoji | `fastqe` | — |
| `fastqscreen` | Screen reads against several genomes for contamination | `fastqscreen_database`, `fastqscreen` | — |
| `fastqutils` | Validate a FASTQ or a pair of mates (fails on malformed input) | `fastqutils_info` | — |
| `fibertoolsrs` | Fiber-seq m6A prediction, nucleosome calling and feature extraction | `fibertoolsrs_predictm6a`, `fibertoolsrs_addnucleosomes`, `fibertoolsrs_extract` | — |
| `filtlong` | Quality filtering of long reads, optionally against Illumina reads | `filtlong` | — |
| `flye` | De novo assembly of PacBio or Nanopore long reads | `flye` | — |
| `fmhfunprofiler` | KEGG ortholog profiling of metagenomes with FracMinHash sketches | `fmhfunprofiler` | — |
| `foldseek` | Protein structure database building and search | `foldseek_createdb`, `foldseek_easysearch` | — |
| `fq` | FASTQ subsampling and validation | `fq_subsample`, `fq_lint` | rnaseq |
| `fqtk` | FASTQ demultiplexing by sample barcodes | `fqtk` | — |
| `gappa` | Analysis of phylogenetic placements: taxonomic assignment, grafted trees, heat trees, merging jplace files | `gappa_examine_assign`, `gappa_examine_graft`, `gappa_examine_heat_tree`, `gappa_edit_merge` | — |
| `gatk4` | GATK: sequence dictionary, interval list and BED conversion, interval-list scattering, MarkDuplicates, base-quality recalibration, HaplotypeCaller, GVCF combining, split-N-CIGAR read splitting, Mutect2 somatic calling and filtering, hard filtering, VCF merging, CNN scoring and tranche filtering | `gatk4_createsequencedictionary`, `gatk4_intervallisttobed`, `gatk4_bedtointervallist`, `gatk4_intervallisttools`, `gatk4_markduplicates`, `gatk4_baserecalibrator`, `gatk4_gatherbqsrreports`, `gatk4_applybqsr`, `gatk4_haplotypecaller`, `gatk4_combinegvcfs`, `gatk4_splitncigarreads`, `gatk4_mutect2`, `gatk4_mergemutectstats`, `gatk4_filtermutectcalls`, `gatk4_variantfiltration`, `gatk4_mergevcfs`, `gatk4_cnnscorevariants`, `gatk4_filtervarianttranches`, `gatk4_preprocessintervals`, `gatk4_annotateintervals`, `gatk4_collectreadcounts`, `gatk4_filterintervals`, `gatk4_createreadcountpanelofnormals`, `gatk4_determinegermlinecontigploidy`, `gatk4_germlinecnvcaller`, `gatk4_genomicsdbimport`, `gatk4_genomicsdbimport_update`, `gatk4_genomicsdbimport_interval_list`, `gatk4_createsomaticpanelofnormals`, `gatk4_indexfeaturefile`, `gatk4_splitintervals` | variant_calling |
| `gawk` | Run an awk program (text or file) over files, decompressing `.gz` inputs and optionally gzipping the output | `gawk` | — |
| `genmap` | Build GenMap's FM index of a FASTA for mappability | `genmap_index` | — |
| `gfaffix` | Collapse redundant shared affixes in a variation graph | `gfaffix` | — |
| `gffread` | Convert GFF to GTF; extract transcript sequences | `gffread_gff_to_gtf`, `gffread_transcripts_fasta` | rnaseq |
| `gmmdemux` | Demultiplex cell-hashing data with Gaussian mixture models (GMM-Demux) | `gmmdemux` | — |
| `gnusort` | Sorting text files (BED, GTF, CSV) with GNU sort | `gnu_sort` | — |
| `graphtyper` | Population-scale small-variant genotyping and region VCF concatenation | `graphtyper_genotype`, `graphtyper_vcfconcatenate` | — |
| `gstama` | TAMA long-read transcript collapse, model merge and FLNC poly(A) cleanup | `gstama_collapse`, `gstama_merge`, `gstama_polyacleanup` | — |
| `gtf2bed` | BED12 gene models from a GTF | `gtf2bed` | rnaseq |
| `gtffilter` | Restrict a GTF to a genome's sequences | `gtf_filter` | rnaseq |
| `gtfgenefilter` | Keep the GTF lines on a genome's sequences (single-cell) | `gtf_gene_filter` | scrnaseq |
| `gunzip` | Decompress a gzipped file | `gunzip` | rnaseq, scrnaseq |
| `gzrt` | Recover the readable data from corrupt gzipped FASTQs | `gzrt` | — |
| `hasheddrops` | Demultiplex cell-hashing data with DropletUtils' hashedDrops (optionally after emptyDrops) | `hasheddrops` | — |
| `hhsuite` | Convert multiple alignments between formats (HH-suite reformat.pl: fas, a2m, a3m, sto, …) | `hhsuite_reformat` | — |
| `hisat2` | Splice-aware short-read alignment: splice sites from a GTF, index build, alignment to BAM | `hisat2_extractsplicesites`, `hisat2_build`, `hisat2_align` | — |
| `hmmer` | Profile HMMs: search a sequence database (hmmsearch), align to a profile (hmmalign), build a profile (hmmbuild), rank hits across profiles (hmmrank), index a sequence file for retrieval (esl-sfetch --index) | `hmmsearch`, `hmmalign`, `hmmbuild`, `hmmrank`, `esl_sfetch_index`, `esl_alimask`, `esl_reformat`, `format_tsv` | — |
| `homer` | HOMER tag directories, peak calling, UCSC tracks, peak-to-BED conversion and peak annotation | `homer_maketagdirectory`, `homer_findpeaks`, `homer_makeucscfile`, `homer_pos2bed`, `homer_annotatepeaks` | — |
| `htodemux` | Demultiplex a Seurat object's hashtag assay with HTODemux | `htodemux` | — |
| `htslib` | bgzip and tabix | `htslib_bgziptabix`, `tabix_tabix`, `tabix_bgzip` | variant_calling |
| `humann3` | Functional profiling of metagenomes with HUMAnN 3: gene families, pathways, regrouping, renormalisation | `humann3_humann`, `humann3_regroup`, `humann3_renorm` | — |
| `ilastik` | Headless pixel classification and multicut segmentation with trained ilastik projects | `ilastik_pixelclassification`, `ilastik_multicut` | — |
| `interproscan` | Protein domain and family annotation | `interproscan` | — |
| `intervals` | Build genomic intervals and split them for scatter-gather | `build_intervals`, `create_intervals_bed` | variant_calling |
| `iqtree` | Maximum-likelihood phylogenies from alignments or trees | `iqtree` | — |
| `isoseq` | Iso-Seq refine: poly(A) trimming and concatemer removal into FLNC reads | `isoseq_refine` | — |
| `kallisto` | Pseudoalignment index and transcript quantification | `kallisto_index`, `kallisto_quant` | — |
| `khmer` | Unique k-mer counting, for effective genome size | `khmer_uniquekmers` | — |
| `kofamscan` | KEGG Orthology assignment with KofamScan's profile HMMs (KEGG's KO data needs a KEGG licence for commercial use) | `kofamscan` | — |
| `kraken2` | Classify reads taxonomically against a Kraken2 database | `kraken2` | — |
| `krona` | Interactive taxonomy charts; build Krona's taxonomy database from NCBI | `krona_ktimporttaxonomy`, `krona_ktupdatetaxonomy` | — |
| `last` | LAST pairwise genome/read alignment: index, train scoring parameters, align, keep best alignments, convert MAF (PSL/SAM/BAM/CRAM/BCF/GFF/…) and draw dot plots | `last_lastdb`, `last_train`, `last_lastal`, `last_split`, `last_mafconvert`, `last_dotplot` | — |
| `leafcutter` | Intron clustering of regtools junctions for differential splicing | `leafcutter_clusterregtools` | — |
| `liftoff` | Lift genome annotations over from a reference assembly to a target | `liftoff` | — |
| `lima` | PacBio primer/barcode removal from CCS reads (BAM, FASTA, FASTQ) | `lima` | — |
| `macs2` | ChIP-seq/ATAC-seq peak calling | `macs2_callpeak` | — |
| `macs3` | Model-based peak calling for ChIP-seq and ATAC-seq | `macs3_callpeak` | — |
| `mafft` | Multiple sequence alignment, optionally adding sequences or fragments | `mafft_align` | — |
| `mageck` | CRISPR screen analysis: sgRNA counting, MLE and RRA gene ranking | `mageck_count`, `mageck_mle`, `mageck_test` | — |
| `mat2json` | Convert MATLAB .mat files to JSON or CSV | `mat2json` | — |
| `mcquant` | Single-cell quantification: per-cell marker intensities and morphology over a segmentation mask | `mcquant` | — |
| `md5sum` | MD5 checksums of files, one per file or combined | `md5sum` | — |
| `medaka` | Consensus polishing of assemblies with nanopore reads | `medaka` | — |
| `megahit` | De novo assembly of short reads, single genomes or metagenomes | `megahit` | — |
| `metaeuk` | Eukaryotic gene prediction in contigs against reference proteins with MetaEuk | `metaeuk_easypredict` | — |
| `metaphlan` | Taxonomic profiling of metagenomes with MetaPhlAn 4 marker genes | `metaphlan_metaphlan` | — |
| `mgikit` | MGI/DNBSEQ FASTQ demultiplexing | `mgikit_demultiplex` | — |
| `mgnifam` | Build protein families (seed and full alignments, HMMs, reference annotations) from a clustering chunk (mgnifam generate_families) | `mgnifam_generatefamilies` | — |
| `mifaser` | Enzyme (EC number) annotation of metagenomic reads with mi-faser | `mifaser` | — |
| `mindagap` | Fill tile grid lines in spatial-transcriptomics panoramas; mark spots duplicated across tile overlaps | `mindagap_mindagap`, `mindagap_duplicatefinder` | — |
| `miniasm` | De novo assembly of long reads from all-vs-all overlaps (no consensus) | `miniasm` | — |
| `minimap2` | Align long or short reads (FASTQ or SAM/BAM/CRAM) to a FASTA or `.mmi` index, writing sorted BAM or PAF; build `.mmi` indexes | `minimap2_index`, `minimap2_align` | — |
| `mirdeep2` | Small-RNA read mapping and known and novel miRNA discovery | `mirdeep2_mapper`, `mirdeep2_mirdeep2` | — |
| `mirtop` | miRNA annotation in mirGFF3, with counts, export and statistics | `mirtop_gff`, `mirtop_counts`, `mirtop_export`, `mirtop_stats` | — |
| `mirtrace` | Quality control for small-RNA sequencing | `mirtrace_qc` | — |
| `misopy` | MISO annotation indexing of a GFF3 (as a .tar.gz) | `misopy_index` | — |
| `mmseqs` | Build, cluster, search, index and tabulate MMseqs2 sequence databases, including taxonomy, expandable profile databases from TSV and public database downloads | `mmseqs_createdb`, `mmseqs_cluster`, `mmseqs_linclust`, `mmseqs_createtsv`, `mmseqs_taxonomy`, `mmseqs_databases`, `mmseqs_createindex`, `mmseqs_tsv2exprofiledb` | — |
| `modkit` | ONT modified-base calling, pileup and BAM tag repair. Oxford Nanopore public licence: in `licensed_env`, not the default `env` | `modkit_callmods`, `modkit_pileup`, `modkit_repair` | — |
| `mosdepth` | Read-depth summaries of BAM/CRAM files | `mosdepth` | variant_calling |
| `msisensorpro` | Scan a reference FASTA for microsatellite sites | `msisensorpro_scan` | — |
| `multiqc` | Aggregate QC reports into one HTML report | `multiqc` | rnaseq, scrnaseq, variant_calling |
| `multiqccustombiotype` | featureCounts biotype counts as MultiQC content | `multiqc_custom_biotype` | rnaseq |
| `multiqcsav` | MultiQC report with Illumina run (SAV/InterOp) metrics | `multiqcsav` | — |
| `multiseqdemux` | Demultiplex a Seurat object's hashtag assay with MULTIseqDemux | `multiseqdemux` | — |
| `nacho` | NanoString nCounter QC reports and normalised counts | `nacho_qc`, `nacho_normalize` | — |
| `nanocomp` | Compare long-read datasets (FASTQ, FASTA, BAM or sequencing summaries): HTML plots and statistics | `nanocomp` | — |
| `nanolyse` | Remove lambda phage reads from Nanopore FASTQ | `nanolyse` | — |
| `nanoplot` | Long-read QC plots and summary statistics from a FASTQ or sequencing summary | `nanoplot` | — |
| `ncbigenomedownload` | Download genome assemblies and annotations from NCBI | `ncbigenomedownload` | — |
| `numorph` | Light-sheet microscopy intensity correction, resampling and 3D U-Net nuclei counting | `numorph_intensity`, `numorph_resample`, `numorph_3dunet` | — |
| `odgi` | Build, sort, transform, lay out, draw and summarise variation graphs | `odgi_build`, `odgi_sort`, `odgi_unchop`, `odgi_squeeze`, `odgi_view`, `odgi_layout`, `odgi_draw`, `odgi_viz`, `odgi_stats` | — |
| `optitype` | HLA class I genotyping from reads aligned to the HLA reference | `optitype` | — |
| `pairix` | Index a bgzipped, sorted contact pairs file (.px2) | `pairix` | — |
| `pairtools` | Parse, sort, deduplicate, flip, restriction-annotate and filter Hi-C contact pairs | `pairtools_parse`, `pairtools_sort`, `pairtools_dedup`, `pairtools_flip`, `pairtools_restrict`, `pairtools_select` | — |
| `pbccs` | HiFi consensus reads (ccs) from one chunk of a PacBio subreads BAM | `pbccs` | — |
| `pbjasmine` | 5mC CpG methylation calling on PacBio HiFi reads | `pbjasmine` | — |
| `pbmm2` | minimap2 alignment for PacBio reads | `pbmm2_align` | — |
| `pcgr` | Download PCGR's reference data bundle (cached, as an archive) | `pcgr_getref` | — |
| `pear` | Paired-end read merging. Non-commercial licence: in `licensed_env`, not the default `env` | `pear` | — |
| `phantompeakqualtools` | ChIP-seq strand cross-correlation quality metrics (NSC, RSC) | `phantompeakqualtools` | — |
| `picard` | Picard MarkDuplicates, CollectMultipleMetrics, MergeSamFiles, AddOrReplaceReadGroups, CreateSequenceDictionary and CollectHsMetrics | `picard_markduplicates`, `picard_collectmultiplemetrics`, `picard_addorreplacereadgroups`, `picard_createsequencedictionary`, `picard_collecthsmetrics`, `picard_mergesamfiles`, `picard_sortsam` | rnaseq |
| `pigz` | Parallel gzip compression and decompression | `pigz_compress`, `pigz_uncompress` | — |
| `pints` | Peak calling for nascent transcription data (PRO-cap, GRO-cap and similar) | `pints_caller` | — |
| `pirate` | Pangenome gene families across prokaryote GFF3 annotations | `pirate` | — |
| `popscle` | Genotype-based demultiplexing: pileup, demuxlet, freemuxlet | `popscle_dscpileup`, `popscle_demuxlet`, `popscle_freemuxlet` | — |
| `porechop` | Adapter trimming of Nanopore reads (Porechop, Porechop_ABI) | `porechop`, `porechop_abi` | — |
| `preseq` | Library complexity from mapped reads: extrapolation and the observed complexity curve | `preseq_lcextrap`, `preseq_ccurve` | — |
| `prodigal` | Prokaryotic gene prediction | `prodigal` | — |
| `prokka` | Annotate a prokaryotic or viral genome assembly | `prokka` | — |
| `pydamage` | Per-contig ancient-DNA damage estimates of an indexed BAM | `pydamage_analyze` | — |
| `qcatch` | Cell calling and cell-level QC of a simpleaf quantification | `qcatch` | scrnaseq |
| `qualimap` | Alignment QC: RNA-seq QC and BAM QC reports | `qualimap_rnaseq`, `qualimap_bamqc` | rnaseq |
| `quarto` | Render a parameterised Quarto, R Markdown or Jupyter notebook to HTML | `quarto_notebook` | — |
| `quast` | Assembly quality assessment, optionally against a reference FASTA and GFF | `quast` | — |
| `racon` | Consensus polishing of long-read assemblies | `racon` | — |
| `rasusa` | Random subsampling of reads to a target depth | `rasusa` | — |
| `raven` | De novo assembly of long uncorrected reads | `raven` | — |
| `rclone` | Copy and verify data between storage remotes (check, checksum, copy) | `rclone_check`, `rclone_checksum`, `rclone_copy` | — |
| `regtools` | Splice junction extraction from RNA-seq BAMs | `regtools_junctionsextract` | — |
| `resolvetaxonomy` | Resolve a taxonomy table from a taxonomy file or sequence headers, and strip headers to bare ids | `resolve_taxonomy` | — |
| `rgi` | Antimicrobial resistance genes from CARD with RGI: database preparation, read mapping, contig/genome prediction (CARD data: no commercial use without a licence) | `rgi_cardannotation`, `rgi_bwt`, `rgi_main` | — |
| `ribodetector` | Remove rRNA reads from FASTQ with a deep-learning classifier (CPU) | `ribodetector` | — |
| `riker` | Alignment QC metrics (alignment, isize, WGS, GC bias, hybcap, error, RNA) in one pass | `riker_multi` | — |
| `rmats` | rMATS' prep step: per-BAM splice graph and junction counts. Non-commercial licence: in `licensed_env`, not the default `env` | `rmats_prep` | — |
| `rsem` | RSEM reference preparation (optionally with a STAR index) and gene/isoform expression estimation | `rsem_preparereference`, `rsem_calculateexpression` | — |
| `rsemmergecounts` | Merge per-sample RSEM results into gene/transcript count and TPM matrices and long tables | `rsem_merge_counts` | — |
| `rseqc` | RNA-seq BAM QC (RSeQC), including transcript integrity (TIN) | `bam_stat`, `infer_experiment`, `inner_distance`, `junction_annotation`, `junction_saturation`, `read_distribution`, `read_duplication`, `tin`, `rseqc` | rnaseq |
| `rtgtools` | Convert FASTA, FASTQ or SAM/BAM into RTG Tools' SDF format | `rtgtools_format` | — |
| `rundbcan` | CAZyme annotation with run_dbcan, and its database download | `rundbcan_database`, `rundbcan_cazymeannotation` | — |
| `rustqc` | One-pass RNA-seq alignment QC (featureCounts, dupRadar, preseq, samtools, RSeQC and Qualimap equivalents) | `rustqc` | — |
| `s4pred` | Protein secondary structure prediction | `s4pred_runmodel` | — |
| `salmon` | Transcript quantification | `salmon_index`, `salmon_quant_reads`, `salmon_quant_bam` | rnaseq |
| `sam2lca` | Lowest-common-ancestor taxonomic assignment of aligned reads | `sam2lca_analyze` | — |
| `samshee` | Illumina sample sheet validation and formatting | `samshee` | — |
| `samtools` | Index, sort, view, merge and summarise alignments; convert them back to FASTQ; convert between BAM and CRAM; build a sequence dictionary; recalculate MD/NM tags, fix mate information and mark duplicates; import FASTQ into unaligned BAM | `samtools_faidx`, `samtools_sort`, `samtools_view`, `samtools_merge`, `samtools_index`, `samtools_stats`, `samtools_flagstat`, `samtools_idxstats`, `samtools_fastq`, `samtools_collatefastq`, `samtools_convert`, `samtools_dict`, `samtools_calmd`, `samtools_fixmate`, `samtools_markdup`, `samtools_import`, `samtools_trimheader` | rnaseq, variant_calling |
| `scanpy` | simpleaf counts to h5ad; concatenate h5ads | `mtx_to_h5ad_simpleaf`, `concat_h5ad`, `scanpy_hashsolo` | scrnaseq |
| `scimap` | Cluster cells of a cell-by-feature table by expression and spatial neighbourhood (AnnData + CSV) | `scimap_mcmicro` | — |
| `seacr` | Peak calling from CUT&RUN or CUT&Tag bedGraphs, against an IgG control or a signal threshold | `seacr_callpeak` | — |
| `seq2hla` | HLA class I/II typing and expression from paired-end RNA-seq reads | `seq2hla` | — |
| `seqcluster` | Collapsing of identical small-RNA reads | `seqcluster_collapse` | — |
| `seqfu` | Sequence-length statistics for FASTA/FASTQ files, as a table and a MultiQC table | `seqfu_stats` | — |
| `seqkit` | FASTA/FASTQ statistics tables; pattern replacement; sequence transformation and filtering; duplicate removal; splitting reads into parts; FASTQ to FASTA conversion; selecting records by name or sequence pattern | `seqkit_stats`, `seqkit_replace`, `seqkit_seq`, `seqkit_rmdup`, `seqkit_split2`, `seqkit_fq2fa`, `seqkit_grep`, `seqkit_head`, `seqkit_sliding` | — |
| `seqtk` | Subsample FASTQ reads (by count or fraction, seeded); quality-trim them; convert or transform FASTA/FASTQ; list the gaps (runs of N) in a FASTA as BED | `seqtk_sample`, `seqtk_trim`, `seqtk_seq`, `seqtk_cutn`, `seqtk_subseq`, `seqtk_mergepe` | — |
| `sequali` | Sequencing-read QC reports (FASTQ or unaligned BAM) | `sequali` | — |
| `seqwish` | Induce a variation graph (GFA) from pairwise alignments | `seqwish_induce` | — |
| `sgdemux` | Singular Genomics G4 FASTQ demultiplexing (G4-only licence) | `sgdemux` | — |
| `simpleaf` | Single-cell quantification with alevin-fry (index, quant) | `simpleaf_index`, `simpleaf_quant` | scrnaseq |
| `smoothxg` | Smooth a variation graph by local partial order realignment | `smoothxg` | — |
| `snpeff` | Variant annotation with snpEff, with a cached download of its database | `snpeff_download`, `snpeff_snpeff` | — |
| `sortmerna` | Filter rRNA reads out of single- or paired-end FASTQ against rRNA databases, or only build their index | `sortmerna` | — |
| `souporcell` | Genotype-free demultiplexing of pooled cells by their variants | `souporcell` | — |
| `sourcepredict` | Predict the source proportions of metagenomic samples from taxon counts | `sourcepredict` | — |
| `sourmash` | MinHash sketches of sequence files, SBT indexes of them, and gather searches | `sourmash_sketch`, `sourmash_index`, `sourmash_gather`, `sourmash_compare` | — |
| `spaceranger` | 10x Space Ranger: build a reference; count Visium and Visium HD spatial gene expression. Licensed: in `licensed_env`, not the default `env` | `spaceranger_mkref`, `spaceranger_count` | — |
| `spades` | De novo genome, RNA, viral and metagenome assembly (Illumina with optional PacBio/Nanopore reads, dataset YAML, custom HMMs) | `spades` | — |
| `sratools` | Download SRA runs (with retries and integrity checks), write the SRA Toolkit settings, dump runs to FASTQ | `sratools_ncbi_settings`, `sratools_prefetch`, `sratools_fasterqdump` | — |
| `star` | Spliced RNA-seq alignment; the STAR genome index version | `star_genome_generate`, `star_align`, `star_indexversion` | rnaseq |
| `stardist` | Nucleus segmentation of 2D images | `stardist` | — |
| `stringtie` | Transcript assembly and quantification; merge assemblies across samples | `stringtie`, `stringtie_merge` | rnaseq |
| `subread` | featureCounts read summarisation | `featurecounts` | rnaseq |
| `summarizedexperiment` | Bundle count matrices into a SummarizedExperiment RDS | `summarized_experiment` | rnaseq |
| `suppa` | Alternative splicing events, PSI per event and isoform, differential splicing and event clustering with SUPPA2 | `suppa_generateevents`, `suppa_psiperevent`, `suppa_psiperisoform`, `suppa_diffsplice`, `suppa_clusterevents` | — |
| `sylph` | Metagenomic profiling by k-mer containment ANI | `sylph_profile` | — |
| `sylphtax` | Taxonomic profiles (MetaPhlAn-style) from sylph results | `sylphtax_taxprof` | — |
| `taxonkit` | Full NCBI-taxonomy lineages of taxon IDs with TaxonKit | `taxonkit_lineage` | — |
| `thermorawfileparser` | Convert Thermo RAW mass-spec files to mzML, MGF or Parquet | `thermorawfileparser` | — |
| `toulligqc` | Nanopore run QC reports from reads, sequencing summaries or BAMs | `toulligqc` | — |
| `transdecoder` | Coding regions in transcripts with TransDecoder | `transdecoder_longorf`, `transdecoder_predict` | — |
| `trimgalore` | Adapter and quality trimming | `trimgalore` | rnaseq |
| `trimmomatic` | Adapter and quality trimming of single- or paired-end Illumina reads | `trimmomatic` | — |
| `trinity` | De novo transcriptome assembly from single-end or paired RNA-seq reads | `trinity` | — |
| `trust4` | TCR/BCR repertoire reconstruction from bulk or single-cell RNA-seq | `trust4` | — |
| `tx2gene` | Transcript-to-gene table from a GTF | `tx2gene` | rnaseq |
| `tximport` | Transcript quantifications to gene count/TPM matrices | `tximport`, `collect_quants` | rnaseq |
| `ucsc` | UCSC utilities: clip bedGraphs, convert bedGraph or wiggle to bigWig, convert GTF to genePred/refFlat | `bedclip`, `bedgraphtobigwig`, `wigtobigwig`, `bigwigaverageoverbed`, `gtftogenepred` | rnaseq |
| `ultra` | uLTRA annotation-guided splice alignment of long transcript reads | `ultra_index`, `ultra_align` | — |
| `umicollapse` | Deduplicate reads by UMI in a sorted BAM or a FASTQ | `umicollapse` | — |
| `umitools` | UMI extraction and deduplication | `umitools_extract`, `umitools_dedup`, `umitools_prepareforrsem` | rnaseq |
| `unicycler` | Hybrid or single-technology bacterial genome assembly | `unicycler` | — |
| `untar` | Extract a tar archive | `untar` | rnaseq |
| `unzip` | Extract a zip (or other 7-Zip-readable) archive into a Dir | `unzip` | — |
| `unzipfiles` | Extract a zip (or other 7-Zip-readable) archive and list every extracted file | `unzipfiles` | — |
| `vcflib` | Filter VCF records/genotypes by field expressions | `vcflib_vcffilter` | — |
| `vcftools` | VCF summaries (Ts/Tv, FILTER, frequencies) | `vcftools` | variant_calling |
| `vireo` | Assign droplets to donors from per-cell allele counts | `vireo` | — |
| `vsearch` | Cluster sequences (fast, size, smallmem, UNOISE); sort sequences by abundance or length | `vsearch_cluster`, `vsearch_sort`, `vsearch_sintax`, `vsearch_usearch_global` | — |
| `wfmash` | All-to-all mapping and alignment of genomes for pangenome graphs | `wfmash` | — |
| `wget` | Download a file over HTTP(S) or FTP with GNU Wget | `wget` | — |
| `wipertools` | Repair malformed FASTQs in parallel chunks: scatter, wipe, gather reads and reports | `wipertools_fastqscatter`, `wipertools_fastqwiper`, `wipertools_fastqgather`, `wipertools_reportgather` | — |
| `xz` | Decompress a single `.xz` file | `xz_decompress` | — |
| `yara` | Read mapping to a reference; mapped reads sorted and indexed | `yara_index`, `yara_mapper` | — |

Each module's docstring describes its tools' inputs, outputs and naming in
detail.
