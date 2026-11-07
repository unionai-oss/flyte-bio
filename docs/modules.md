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
| `bakta` | Annotate a bacterial genome, plasmid or MAG; download its database | `bakta_db_download`, `bakta` | — |
| `bamtools` | BAM statistics; split BAMs by reference, mapping, pairing or tag, as a `.tar.gz` of BAMs | `bamtools_stats`, `bamtools_split` | — |
| `bbmap` | BBTools: BBSplit binning of reads by reference; BBDuk trimming; BBNorm read-depth normalisation; FilterByName read selection; Repair re-pairing of out-of-order mates | `bbsplit_index`, `bbsplit`, `bbduk`, `bbnorm`, `bbmap_filterbyname`, `bbmap_repair` | rnaseq |
| `bcftools` | VCF/BCF utilities: statistics, sort, normalise, view, filter, annotate, merge, mpileup calling, index, concatenate, reheader, consensus sequence, set genotypes | `bcftools_stats`, `bcftools_sort`, `bcftools_norm`, `bcftools_view`, `bcftools_filter`, `bcftools_annotate`, `bcftools_merge`, `bcftools_mpileup`, `bcftools_index`, `bcftools_concat`, `bcftools_reheader`, `bcftools_consensus`, `bcftools_pluginsetgt` | variant_calling |
| `bedops` | GTF to BED conversion | `bedops_gtf2bed` | — |
| `bedtools` | Genome arithmetic; extract or mask a BED's intervals in a FASTA; convert a BAM to BED; complement a BED against a genome | `bedtools_intersect`, `bedtools_sort`, `bedtools_merge`, `bedtools_genomecov`, `bedtools_getfasta`, `bedtools_maskfasta`, `bedtools_bamtobed`, `bedtools_complement` | rnaseq |
| `blast` | Build BLAST databases and search nucleotide sequences with blastn | `blast_makeblastdb`, `blast_blastn` | — |
| `bowtie` | Short-read alignment with Bowtie 1: index build, alignment to BAM | `bowtie1_build`, `bowtie1_align` | — |
| `bowtie2` | Bowtie 2 short-read alignment: build an index, align reads through samtools to BAM/SAM/CRAM, with the log and optional unaligned reads | `bowtie2_build`, `bowtie2_align` | — |
| `bracken` | Abundance re-estimation from a Kraken2 report | `bracken` | — |
| `busco` | Assembly, proteome and transcriptome completeness against a BUSCO lineage; download lineages | `busco_download`, `busco` | — |
| `bwa` | BWA short-read alignment: index, BWA-MEM, and BWA-backtrack (`aln`, then `samse`/`sampe` to sorted BAM) | `bwa_index`, `bwa_mem`, `bwa_aln`, `bwa_samse`, `bwa_sampe` | variant_calling |
| `bwamem2` | BWA-MEM2 index and short-read alignment, piped through samtools sort or view | `bwamem2_index`, `bwamem2_mem` | — |
| `cat` | Concatenate a sample's FASTQs; concatenate any files, (de)compressing by extension | `cat_fastq`, `cat_cat`, `concatenate` | rnaseq |
| `catadditionalfasta` | Append extra sequences (e.g. spike-ins) to a genome and GTF | `cat_additional_fasta` | rnaseq |
| `cellbender` | Remove ambient RNA and empty droplets (CPU or GPU) | `cellbender_removebackground` | scrnaseq |
| `cellpose` | Cell segmentation of microscopy images | `cellpose` | — |
| `cellranger` | 10x Cell Ranger V(D)J reference building and immune profiling. Licensed: in `licensed_env`, not the default `env` | `cellranger_mkvdjref`, `cellranger_vdj` | — |
| `chromsizes` | Chromosome sizes and FASTA index (`.fai`, `.gzi` if bgzipped) | `chrom_sizes` | — |
| `clustalo` | Multiple sequence alignment and guide trees with Clustal Omega | `clustalo_align`, `clustalo_guidetree` | — |
| `cooler` | Hi-C contact matrices: bin pairs into a .cool, make genomic bins, balance, dump tables, zoomify to a multi-resolution .mcool | `cooler_cload`, `cooler_makebins`, `cooler_balance`, `cooler_dump`, `cooler_zoomify` | — |
| `csvtk` | Join tables on key fields; concatenate tables, converting the delimiter | `csvtk_join`, `csvtk_concat` | — |
| `cutadapt` | Adapter, primer and quality trimming of FASTQ reads (single- or paired-end) | `cutadapt` | — |
| `damageprofiler` | Ancient-DNA damage patterns (misincorporation, read lengths) of a BAM | `damageprofiler` | — |
| `deepcell` | Whole-cell and nuclear segmentation with Mesmer. Non-commercial licence: in `licensed_env`, not the default `env` | `deepcell_mesmer` | — |
| `deeptools` | Coverage matrices, heatmaps and fingerprint plots for ChIP/ATAC-seq; BAM bin summaries with correlation and PCA plots | `deeptools_computematrix`, `deeptools_plotheatmap`, `deeptools_plotfingerprint`, `deeptools_multibamsummary`, `deeptools_plotcorrelation`, `deeptools_plotpca` | — |
| `deepvariant` | Germline small-variant calling with DeepVariant (VCF and gVCF) | `deepvariant_rundeepvariant` | — |
| `deseq2_qc` | Sample-level QC of a count matrix with DESeq2 (PCA, distances) | `deseq2_qc` | rnaseq |
| `diamond` | Cluster protein sequences (FASTA or DIAMOND database) | `diamond_cluster` | — |
| `dragmap` | Short-read alignment with DRAGMAP: hash table build, alignment to BAM/CRAM | `dragmap_hashtable`, `dragmap_align` | — |
| `duckdb` | Convert a CSV/TSV (optionally compressed) to Parquet | `duckdb_table2parquet` | — |
| `dupradar` | Duplication rate vs. expression for RNA-seq BAMs | `dupradar` | rnaseq |
| `ensemblvep` | Variant annotation with Ensembl VEP, with its cache download | `ensemblvep_download`, `ensemblvep_vep` | — |
| `evigene` | Reduce an over-assembled transcriptome to non-redundant transcripts, CDS and proteins | `evigene_tr2aacds` | — |
| `fastp` | Adapter and quality trimming with QC reports (single-end, paired-end, interleaved) | `fastp` | — |
| `fastqc` | Read quality reports | `fastqc`, `fastqc_reads` | rnaseq, scrnaseq, variant_calling |
| `fastqdl` | Download FASTQ files and run metadata from ENA or SRA | `fastqdl` | — |
| `fastqutils` | Validate a FASTQ or a pair of mates (fails on malformed input) | `fastqutils_info` | — |
| `filtlong` | Quality filtering of long reads, optionally against Illumina reads | `filtlong` | — |
| `flye` | De novo assembly of PacBio or Nanopore long reads | `flye` | — |
| `foldseek` | Protein structure database building and search | `foldseek_createdb`, `foldseek_easysearch` | — |
| `fq` | FASTQ subsampling and validation | `fq_subsample`, `fq_lint` | rnaseq |
| `gatk4` | GATK: sequence dictionary, interval list and BED conversion, interval-list scattering, MarkDuplicates, base-quality recalibration, HaplotypeCaller, GVCF combining, split-N-CIGAR read splitting, Mutect2 somatic calling and filtering, hard filtering, VCF merging, CNN scoring and tranche filtering | `gatk4_createsequencedictionary`, `gatk4_intervallisttobed`, `gatk4_bedtointervallist`, `gatk4_intervallisttools`, `gatk4_markduplicates`, `gatk4_baserecalibrator`, `gatk4_gatherbqsrreports`, `gatk4_applybqsr`, `gatk4_haplotypecaller`, `gatk4_combinegvcfs`, `gatk4_splitncigarreads`, `gatk4_mutect2`, `gatk4_mergemutectstats`, `gatk4_filtermutectcalls`, `gatk4_variantfiltration`, `gatk4_mergevcfs`, `gatk4_cnnscorevariants`, `gatk4_filtervarianttranches` | variant_calling |
| `gawk` | Run an awk program (text or file) over files, decompressing `.gz` inputs and optionally gzipping the output | `gawk` | — |
| `gffread` | Convert GFF to GTF; extract transcript sequences | `gffread_gff_to_gtf`, `gffread_transcripts_fasta` | rnaseq |
| `gtf2bed` | BED12 gene models from a GTF | `gtf2bed` | rnaseq |
| `gtffilter` | Restrict a GTF to a genome's sequences | `gtf_filter` | rnaseq |
| `gtfgenefilter` | Keep the GTF lines on a genome's sequences (single-cell) | `gtf_gene_filter` | scrnaseq |
| `gunzip` | Decompress a gzipped file | `gunzip` | rnaseq, scrnaseq |
| `gzrt` | Recover the readable data from corrupt gzipped FASTQs | `gzrt` | — |
| `hisat2` | Splice-aware short-read alignment: splice sites from a GTF, index build, alignment to BAM | `hisat2_extractsplicesites`, `hisat2_build`, `hisat2_align` | — |
| `hmmer` | Profile HMMs: search a sequence database (hmmsearch), align to a profile (hmmalign), build a profile (hmmbuild), rank hits across profiles (hmmrank) | `hmmsearch`, `hmmalign`, `hmmbuild`, `hmmrank` | — |
| `htslib` | bgzip and tabix | `htslib_bgziptabix`, `tabix_tabix`, `tabix_bgzip` | variant_calling |
| `ilastik` | Headless pixel classification and multicut segmentation with trained ilastik projects | `ilastik_pixelclassification`, `ilastik_multicut` | — |
| `interproscan` | Protein domain and family annotation | `interproscan` | — |
| `intervals` | Build genomic intervals and split them for scatter-gather | `build_intervals`, `create_intervals_bed` | variant_calling |
| `kallisto` | Build a kallisto pseudoalignment index from a transcript FASTA | `kallisto_index` | — |
| `kraken2` | Classify reads taxonomically against a Kraken2 database | `kraken2` | — |
| `krona` | Interactive taxonomy charts; build Krona's taxonomy database from NCBI | `krona_ktimporttaxonomy`, `krona_ktupdatetaxonomy` | — |
| `macs2` | ChIP-seq/ATAC-seq peak calling | `macs2_callpeak` | — |
| `mafft` | Multiple sequence alignment, optionally adding sequences or fragments | `mafft_align` | — |
| `mat2json` | Convert MATLAB .mat files to JSON or CSV | `mat2json` | — |
| `md5sum` | MD5 checksums of files, one per file or combined | `md5sum` | — |
| `megahit` | De novo assembly of short reads, single genomes or metagenomes | `megahit` | — |
| `mindagap` | Fill tile grid lines in spatial-transcriptomics panoramas; mark spots duplicated across tile overlaps | `mindagap_mindagap`, `mindagap_duplicatefinder` | — |
| `minimap2` | Align long or short reads (FASTQ or SAM/BAM/CRAM) to a FASTA or `.mmi` index, writing sorted BAM or PAF; build `.mmi` indexes | `minimap2_index`, `minimap2_align` | — |
| `mmseqs` | Build, cluster, search, index and tabulate MMseqs2 sequence databases, including taxonomy, expandable profile databases from TSV and public database downloads | `mmseqs_createdb`, `mmseqs_cluster`, `mmseqs_linclust`, `mmseqs_createtsv`, `mmseqs_taxonomy`, `mmseqs_databases`, `mmseqs_createindex`, `mmseqs_tsv2exprofiledb` | — |
| `mosdepth` | Read-depth summaries of BAM/CRAM files | `mosdepth` | variant_calling |
| `msisensorpro` | Scan a reference FASTA for microsatellite sites | `msisensorpro_scan` | — |
| `multiqc` | Aggregate QC reports into one HTML report | `multiqc` | rnaseq, scrnaseq, variant_calling |
| `multiqccustombiotype` | featureCounts biotype counts as MultiQC content | `multiqc_custom_biotype` | rnaseq |
| `nacho` | NanoString nCounter QC reports and normalised counts | `nacho_qc`, `nacho_normalize` | — |
| `nanocomp` | Compare long-read datasets (FASTQ, FASTA, BAM or sequencing summaries): HTML plots and statistics | `nanocomp` | — |
| `nanolyse` | Remove lambda phage reads from Nanopore FASTQ | `nanolyse` | — |
| `nanoplot` | Long-read QC plots and summary statistics from a FASTQ or sequencing summary | `nanoplot` | — |
| `ncbigenomedownload` | Download genome assemblies and annotations from NCBI | `ncbigenomedownload` | — |
| `numorph` | Light-sheet microscopy intensity correction, resampling and 3D U-Net nuclei counting | `numorph_intensity`, `numorph_resample`, `numorph_3dunet` | — |
| `optitype` | HLA class I genotyping from reads aligned to the HLA reference | `optitype` | — |
| `pcgr` | Download PCGR's reference data bundle (cached, as an archive) | `pcgr_getref` | — |
| `picard` | Picard MarkDuplicates, CollectMultipleMetrics, AddOrReplaceReadGroups, CreateSequenceDictionary and CollectHsMetrics | `picard_markduplicates`, `picard_collectmultiplemetrics`, `picard_addorreplacereadgroups`, `picard_createsequencedictionary`, `picard_collecthsmetrics` | rnaseq |
| `pigz` | Parallel gzip compression and decompression | `pigz_compress`, `pigz_uncompress` | — |
| `porechop` | Adapter trimming of Nanopore reads (Porechop, Porechop_ABI) | `porechop`, `porechop_abi` | — |
| `preseq` | Library complexity from mapped reads: extrapolation and the observed complexity curve | `preseq_lcextrap`, `preseq_ccurve` | — |
| `prodigal` | Prokaryotic gene prediction | `prodigal` | — |
| `prokka` | Annotate a prokaryotic or viral genome assembly | `prokka` | — |
| `pydamage` | Per-contig ancient-DNA damage estimates of an indexed BAM | `pydamage_analyze` | — |
| `qcatch` | Cell calling and cell-level QC of a simpleaf quantification | `qcatch` | scrnaseq |
| `qualimap` | Alignment QC: RNA-seq QC and BAM QC reports | `qualimap_rnaseq`, `qualimap_bamqc` | rnaseq |
| `quarto` | Render a parameterised Quarto, R Markdown or Jupyter notebook to HTML | `quarto_notebook` | — |
| `quast` | Assembly quality assessment, optionally against a reference FASTA and GFF | `quast` | — |
| `rclone` | Copy and verify data between storage remotes (check, checksum, copy) | `rclone_check`, `rclone_checksum`, `rclone_copy` | — |
| `rsem` | Prepare an RSEM reference (optionally with a STAR index) from a genome FASTA and GTF | `rsem_preparereference` | — |
| `rseqc` | RNA-seq BAM QC (RSeQC), including transcript integrity (TIN) | `bam_stat`, `infer_experiment`, `inner_distance`, `junction_annotation`, `junction_saturation`, `read_distribution`, `read_duplication`, `tin`, `rseqc` | rnaseq |
| `rtgtools` | Convert FASTA, FASTQ or SAM/BAM into RTG Tools' SDF format | `rtgtools_format` | — |
| `s4pred` | Protein secondary structure prediction | `s4pred_runmodel` | — |
| `salmon` | Transcript quantification | `salmon_index`, `salmon_quant_reads`, `salmon_quant_bam` | rnaseq |
| `sam2lca` | Lowest-common-ancestor taxonomic assignment of aligned reads | `sam2lca_analyze` | — |
| `samtools` | Index, sort, view, merge and summarise alignments; convert them back to FASTQ; convert between BAM and CRAM; build a sequence dictionary; recalculate MD/NM tags, fix mate information and mark duplicates | `samtools_faidx`, `samtools_sort`, `samtools_view`, `samtools_merge`, `samtools_index`, `samtools_stats`, `samtools_flagstat`, `samtools_idxstats`, `samtools_fastq`, `samtools_collatefastq`, `samtools_convert`, `samtools_dict`, `samtools_calmd`, `samtools_fixmate`, `samtools_markdup` | rnaseq, variant_calling |
| `scanpy` | simpleaf counts to h5ad; concatenate h5ads | `mtx_to_h5ad_simpleaf`, `concat_h5ad` | scrnaseq |
| `seacr` | Peak calling from CUT&RUN or CUT&Tag bedGraphs, against an IgG control or a signal threshold | `seacr_callpeak` | — |
| `seq2hla` | HLA class I/II typing and expression from paired-end RNA-seq reads | `seq2hla` | — |
| `seqfu` | Sequence-length statistics for FASTA/FASTQ files, as a table and a MultiQC table | `seqfu_stats` | — |
| `seqkit` | FASTA/FASTQ statistics tables; pattern replacement; sequence transformation and filtering; duplicate removal; splitting reads into parts | `seqkit_stats`, `seqkit_replace`, `seqkit_seq`, `seqkit_rmdup`, `seqkit_split2` | — |
| `seqtk` | Subsample FASTQ reads (by count or fraction, seeded); quality-trim them | `seqtk_sample`, `seqtk_trim` | — |
| `simpleaf` | Single-cell quantification with alevin-fry (index, quant) | `simpleaf_index`, `simpleaf_quant` | scrnaseq |
| `snpeff` | Variant annotation with snpEff, with a cached download of its database | `snpeff_download`, `snpeff_snpeff` | — |
| `sortmerna` | Filter rRNA reads out of single- or paired-end FASTQ against rRNA databases, or only build their index | `sortmerna` | — |
| `sourcepredict` | Predict the source proportions of metagenomic samples from taxon counts | `sourcepredict` | — |
| `sourmash` | MinHash sketches of sequence files | `sourmash_sketch` | — |
| `spaceranger` | 10x Space Ranger: build a reference; count Visium and Visium HD spatial gene expression. Licensed: in `licensed_env`, not the default `env` | `spaceranger_mkref`, `spaceranger_count` | — |
| `spades` | De novo genome, RNA, viral and metagenome assembly (Illumina with optional PacBio/Nanopore reads, dataset YAML, custom HMMs) | `spades` | — |
| `sratools` | Download SRA runs (with retries and integrity checks), write the SRA Toolkit settings, dump runs to FASTQ | `sratools_ncbi_settings`, `sratools_prefetch`, `sratools_fasterqdump` | — |
| `star` | Spliced RNA-seq alignment; the STAR genome index version | `star_genome_generate`, `star_align`, `star_indexversion` | rnaseq |
| `stardist` | Nucleus segmentation of 2D images | `stardist` | — |
| `stringtie` | Transcript assembly and quantification; merge assemblies across samples | `stringtie`, `stringtie_merge` | rnaseq |
| `subread` | featureCounts read summarisation | `featurecounts` | rnaseq |
| `summarizedexperiment` | Bundle count matrices into a SummarizedExperiment RDS | `summarized_experiment` | rnaseq |
| `thermorawfileparser` | Convert Thermo RAW mass-spec files to mzML, MGF or Parquet | `thermorawfileparser` | — |
| `toulligqc` | Nanopore run QC reports from reads, sequencing summaries or BAMs | `toulligqc` | — |
| `trimgalore` | Adapter and quality trimming | `trimgalore` | rnaseq |
| `trimmomatic` | Adapter and quality trimming of single- or paired-end Illumina reads | `trimmomatic` | — |
| `trinity` | De novo transcriptome assembly from single-end or paired RNA-seq reads | `trinity` | — |
| `trust4` | TCR/BCR repertoire reconstruction from bulk or single-cell RNA-seq | `trust4` | — |
| `tx2gene` | Transcript-to-gene table from a GTF | `tx2gene` | rnaseq |
| `tximport` | Transcript quantifications to gene count/TPM matrices | `tximport`, `collect_quants` | rnaseq |
| `ucsc` | UCSC utilities: clip bedGraphs, convert bedGraph or wiggle to bigWig, convert GTF to genePred/refFlat | `bedclip`, `bedgraphtobigwig`, `wigtobigwig`, `gtftogenepred` | rnaseq |
| `umicollapse` | Deduplicate reads by UMI in a sorted BAM or a FASTQ | `umicollapse` | — |
| `umitools` | UMI extraction and deduplication | `umitools_extract`, `umitools_dedup`, `umitools_prepareforrsem` | rnaseq |
| `untar` | Extract a tar archive | `untar` | rnaseq |
| `unzip` | Extract a zip (or other 7-Zip-readable) archive into a Dir | `unzip` | — |
| `unzipfiles` | Extract a zip (or other 7-Zip-readable) archive and list every extracted file | `unzipfiles` | — |
| `vcftools` | VCF summaries (Ts/Tv, FILTER, frequencies) | `vcftools` | variant_calling |
| `vsearch` | Cluster sequences (fast, size, smallmem, UNOISE) | `vsearch_cluster` | — |
| `wipertools` | Repair malformed FASTQs in parallel chunks: scatter, wipe, gather reads and reports | `wipertools_fastqscatter`, `wipertools_fastqwiper`, `wipertools_fastqgather`, `wipertools_reportgather` | — |
| `xz` | Decompress a single `.xz` file | `xz_decompress` | — |
| `yara` | Read mapping to a reference; mapped reads sorted and indexed | `yara_index`, `yara_mapper` | — |

Each module's docstring describes its tools' inputs, outputs and naming in
detail.
