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
| `anndata` | Keep an h5ad's cells with listed barcodes | `anndata_barcodes` | scrnaseq |
| `anndatar` | h5ad to Seurat and SingleCellExperiment RDS | `anndatar_convert` | scrnaseq |
| `bbmap` | BBTools: BBSplit binning of reads by the reference they map to best; BBDuk quality, adapter and contaminant trimming | `bbsplit_index`, `bbsplit`, `bbduk` | rnaseq |
| `bcftools` | VCF/BCF utilities: statistics, sort, normalise, view, filter, annotate, merge, mpileup calling | `bcftools_stats`, `bcftools_sort`, `bcftools_norm`, `bcftools_view`, `bcftools_filter`, `bcftools_annotate`, `bcftools_merge`, `bcftools_mpileup` | variant_calling |
| `bedtools` | Genome arithmetic | `bedtools_intersect`, `bedtools_sort`, `bedtools_merge`, `bedtools_genomecov` | rnaseq |
| `bowtie2` | Bowtie 2 short-read alignment: build an index, align reads through samtools to BAM/SAM/CRAM, with the log and optional unaligned reads | `bowtie2_build`, `bowtie2_align` | — |
| `bwa` | BWA short-read alignment | `bwa_index`, `bwa_mem` | variant_calling |
| `bwamem2` | BWA-MEM2 index and short-read alignment, piped through samtools sort or view | `bwamem2_index`, `bwamem2_mem` | — |
| `cat` | Concatenate a sample's FASTQs; concatenate any files, (de)compressing by extension | `cat_fastq`, `cat_cat`, `concatenate` | rnaseq |
| `catadditionalfasta` | Append extra sequences (e.g. spike-ins) to a genome and GTF | `cat_additional_fasta` | rnaseq |
| `cellbender` | Remove ambient RNA and empty droplets (CPU or GPU) | `cellbender_removebackground` | scrnaseq |
| `chromsizes` | Chromosome sizes and FASTA index (`.fai`, `.gzi` if bgzipped) | `chrom_sizes` | — |
| `cutadapt` | Adapter, primer and quality trimming of FASTQ reads (single- or paired-end) | `cutadapt` | — |
| `deseq2_qc` | Sample-level QC of a count matrix with DESeq2 (PCA, distances) | `deseq2_qc` | rnaseq |
| `dupradar` | Duplication rate vs. expression for RNA-seq BAMs | `dupradar` | rnaseq |
| `fastp` | Adapter and quality trimming with QC reports (single-end, paired-end, interleaved) | `fastp` | — |
| `fastqc` | Read quality reports | `fastqc`, `fastqc_reads` | rnaseq, scrnaseq, variant_calling |
| `fq` | FASTQ subsampling and validation | `fq_subsample`, `fq_lint` | rnaseq |
| `gatk4` | GATK: sequence dictionary, interval list and BED conversion, MarkDuplicates, base-quality recalibration, HaplotypeCaller, VCF merging, CNN scoring and tranche filtering | `gatk4_createsequencedictionary`, `gatk4_intervallisttobed`, `gatk4_bedtointervallist`, `gatk4_markduplicates`, `gatk4_baserecalibrator`, `gatk4_gatherbqsrreports`, `gatk4_applybqsr`, `gatk4_haplotypecaller`, `gatk4_mergevcfs`, `gatk4_cnnscorevariants`, `gatk4_filtervarianttranches` | variant_calling |
| `gawk` | Run an awk program (text or file) over files, decompressing `.gz` inputs and optionally gzipping the output | `gawk` | — |
| `gffread` | Convert GFF to GTF; extract transcript sequences | `gffread_gff_to_gtf`, `gffread_transcripts_fasta` | rnaseq |
| `gtf2bed` | BED12 gene models from a GTF | `gtf2bed` | rnaseq |
| `gtffilter` | Restrict a GTF to a genome's sequences | `gtf_filter` | rnaseq |
| `gtfgenefilter` | Keep the GTF lines on a genome's sequences (single-cell) | `gtf_gene_filter` | scrnaseq |
| `gunzip` | Decompress a gzipped file | `gunzip` | rnaseq, scrnaseq |
| `htslib` | bgzip and tabix | `htslib_bgziptabix`, `tabix_tabix`, `tabix_bgzip` | variant_calling |
| `intervals` | Build genomic intervals and split them for scatter-gather | `build_intervals`, `create_intervals_bed` | variant_calling |
| `kraken2` | Classify reads taxonomically against a Kraken2 database | `kraken2` | — |
| `minimap2` | Align long or short reads (FASTQ or SAM/BAM/CRAM) to a FASTA or `.mmi` index, writing sorted BAM or PAF; build `.mmi` indexes | `minimap2_index`, `minimap2_align` | — |
| `mosdepth` | Read-depth summaries of BAM/CRAM files | `mosdepth` | variant_calling |
| `multiqc` | Aggregate QC reports into one HTML report | `multiqc` | rnaseq, scrnaseq, variant_calling |
| `multiqccustombiotype` | featureCounts biotype counts as MultiQC content | `multiqc_custom_biotype` | rnaseq |
| `nanoplot` | Long-read QC plots and summary statistics from a FASTQ or sequencing summary | `nanoplot` | — |
| `picard` | Picard MarkDuplicates and CollectMultipleMetrics | `picard_markduplicates`, `picard_collectmultiplemetrics` | rnaseq |
| `pigz` | Parallel gzip compression and decompression | `pigz_compress`, `pigz_uncompress` | — |
| `preseq` | Library complexity extrapolation (`preseq lc_extrap`) from mapped reads | `preseq_lcextrap` | — |
| `qcatch` | Cell calling and cell-level QC of a simpleaf quantification | `qcatch` | scrnaseq |
| `qualimap` | RNA-seq alignment QC | `qualimap_rnaseq` | rnaseq |
| `rseqc` | RNA-seq BAM QC (RSeQC) | `bam_stat`, `infer_experiment`, `inner_distance`, `junction_annotation`, `junction_saturation`, `read_distribution`, `read_duplication`, `rseqc` | rnaseq |
| `salmon` | Transcript quantification | `salmon_index`, `salmon_quant_reads`, `salmon_quant_bam` | rnaseq |
| `samtools` | Index, sort, view, merge and summarise alignments; convert them back to FASTQ | `samtools_faidx`, `samtools_sort`, `samtools_view`, `samtools_merge`, `samtools_index`, `samtools_stats`, `samtools_flagstat`, `samtools_idxstats`, `samtools_fastq`, `samtools_collatefastq` | rnaseq, variant_calling |
| `scanpy` | simpleaf counts to h5ad; concatenate h5ads | `mtx_to_h5ad_simpleaf`, `concat_h5ad` | scrnaseq |
| `seqkit` | FASTA/FASTQ statistics tables; pattern replacement in names or sequences | `seqkit_stats`, `seqkit_replace` | — |
| `simpleaf` | Single-cell quantification with alevin-fry (index, quant) | `simpleaf_index`, `simpleaf_quant` | scrnaseq |
| `star` | Spliced RNA-seq alignment | `star_genome_generate`, `star_align` | rnaseq |
| `stringtie` | Transcript assembly and quantification | `stringtie` | rnaseq |
| `subread` | featureCounts read summarisation | `featurecounts` | rnaseq |
| `summarizedexperiment` | Bundle count matrices into a SummarizedExperiment RDS | `summarized_experiment` | rnaseq |
| `trimgalore` | Adapter and quality trimming | `trimgalore` | rnaseq |
| `tx2gene` | Transcript-to-gene table from a GTF | `tx2gene` | rnaseq |
| `tximport` | Transcript quantifications to gene count/TPM matrices | `tximport`, `collect_quants` | rnaseq |
| `ucsc` | UCSC utilities: clip bedGraphs, convert to bigWig | `bedclip`, `bedgraphtobigwig` | rnaseq |
| `umitools` | UMI extraction and deduplication | `umitools_extract`, `umitools_dedup`, `umitools_prepareforrsem` | rnaseq |
| `untar` | Extract a tar archive | `untar` | rnaseq |
| `unzip` | Extract a zip (or other 7-Zip-readable) archive into a Dir | `unzip` | — |
| `vcftools` | VCF summaries (Ts/Tv, FILTER, frequencies) | `vcftools` | variant_calling |

Each module's docstring describes its tools' inputs, outputs and naming in
detail.
