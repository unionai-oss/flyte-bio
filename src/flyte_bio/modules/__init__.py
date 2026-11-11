"""flyte_bio.modules — typed shell-task wrappers around bio CLI tools.

Each submodule wraps one tool family, sharing a single biocontainer image
across its tasks:

- :mod:`flyte_bio.modules.aardvark` — aardvark small-variant benchmarking.
- :mod:`flyte_bio.modules.abacas` — ABACAS contig ordering.
- :mod:`flyte_bio.modules.adapterremoval` — AdapterRemoval trimming.
- :mod:`flyte_bio.modules.affy` — affy RMA normalisation.
- :mod:`flyte_bio.modules.agat` — AGAT annotation conversion.
- :mod:`flyte_bio.modules.amulety` — AIRR translation and sequence embeddings.
- :mod:`flyte_bio.modules.anndata` — subset an h5ad to listed barcodes.
- :mod:`flyte_bio.modules.anndatar` — h5ad to Seurat / SingleCellExperiment RDS.
- :mod:`flyte_bio.modules.annotsv` — AnnotSV SV annotation (annotation data licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.aria2` — file downloads.
- :mod:`flyte_bio.modules.arriba` — Arriba fusion calling.
- :mod:`flyte_bio.modules.art` — Illumina read simulation.
- :mod:`flyte_bio.modules.artic` — ARTIC amplicon pipeline.
- :mod:`flyte_bio.modules.ashlar` — ASHLAR image stitching and registration.
- :mod:`flyte_bio.modules.assemblyscan` — assembly summary statistics.
- :mod:`flyte_bio.modules.ataqv` — ATAC-seq quality metrics and interactive reports.
- :mod:`flyte_bio.modules.atlasgeneannotationmanipulation` — GTF feature annotation tables.
- :mod:`flyte_bio.modules.autocycler` — Autocycler consensus assembly.
- :mod:`flyte_bio.modules.backsub` — multiplexed image background subtraction.
- :mod:`flyte_bio.modules.bakta` — bacterial genome annotation (database download, annotate).
- :mod:`flyte_bio.modules.bamtools` — BAM statistics and splitting.
- :mod:`flyte_bio.modules.bandage` — Bandage graph images.
- :mod:`flyte_bio.modules.barrnap` — barrnap rRNA prediction.
- :mod:`flyte_bio.modules.bases2fastq` — Element run-folder demultiplexing (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.basicpy` — BaSiC illumination profiles.
- :mod:`flyte_bio.modules.bbmap` — BBSplit (index build + read binning).
- :mod:`flyte_bio.modules.bcftools` — VCF/BCF utilities (stats, sort, norm, view, filter, annotate, merge, mpileup).
- :mod:`flyte_bio.modules.bcl2fastq` — Illumina bcl2fastq (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.bclconvert` — Illumina BCL Convert (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.beagle5` — Beagle 5 phasing and imputation.
- :mod:`flyte_bio.modules.bedops` — GTF to BED conversion.
- :mod:`flyte_bio.modules.bedtools` — genome arithmetic (intersect, sort,
  merge, genomecov).
- :mod:`flyte_bio.modules.bff` — cellhashR BFF cell-hashing demultiplexing.
- :mod:`flyte_bio.modules.bftools` — Bio-Formats tools.
- :mod:`flyte_bio.modules.bioawk` — awk with built-in parsing of sequence and interval formats.
- :mod:`flyte_bio.modules.bismark` — Bismark bisulfite alignment and methylation calls.
- :mod:`flyte_bio.modules.blast` — BLAST databases and blastn search.
- :mod:`flyte_bio.modules.blat` — BLAT alignment (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.bowtie` — Bowtie 1 short-read alignment (build, align).
- :mod:`flyte_bio.modules.bowtie2` — Bowtie 2 short-read alignment (build, align).
- :mod:`flyte_bio.modules.bracken` — abundance re-estimation from Kraken2 reports.
- :mod:`flyte_bio.modules.busco` — assembly and proteome completeness (lineage download, assess).
- :mod:`flyte_bio.modules.bwa` — BWA short-read alignment (index).
- :mod:`flyte_bio.modules.bwamem2` — BWA-MEM2 short-read alignment (index, mem).
- :mod:`flyte_bio.modules.bwameth` — bwa-meth bisulfite alignment.
- :mod:`flyte_bio.modules.canu` — Canu long-read assembly.
- :mod:`flyte_bio.modules.carveme` — CarveMe metabolic models.
- :mod:`flyte_bio.modules.cellranger` — 10x Cell Ranger V(D)J (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.cat` — file concatenation (``cat_fastq``, ``cat_cat``, ``concatenate``).
- :mod:`flyte_bio.modules.deepcell` — Mesmer cell segmentation (non-commercial; in ``licensed_env``).
- :mod:`flyte_bio.modules.deseq2_qc` — DESeq2 PCA / sample-distance QC.
- :mod:`flyte_bio.modules.dupradar` — duplication rate vs expression QC.
- :mod:`flyte_bio.modules.fastp` — read trimming with QC reports.
- :mod:`flyte_bio.modules.fastqc` — read quality-control reports.
- :mod:`flyte_bio.modules.catadditionalfasta` — append an extra FASTA + GTF.
- :mod:`flyte_bio.modules.catpack` — CAT/BAT contig and bin classification.
- :mod:`flyte_bio.modules.cdhit` — CD-HIT clustering.
- :mod:`flyte_bio.modules.cellbender` — ambient RNA and empty droplet removal (CPU or GPU).
- :mod:`flyte_bio.modules.cellpose` — cell segmentation of microscopy images.
- :mod:`flyte_bio.modules.cellrangerarc` — 10x Cell Ranger ARC (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.cellsnp` — cellsnp-lite per-cell allele counts.
- :mod:`flyte_bio.modules.centrifuge` — Centrifuge index building.
- :mod:`flyte_bio.modules.centrifuger` — Centrifuger index building.
- :mod:`flyte_bio.modules.checkm2` — CheckM2 genome completeness.
- :mod:`flyte_bio.modules.checkqc` — Illumina run-folder QC thresholds.
- :mod:`flyte_bio.modules.checkv` — CheckV viral genome quality.
- :mod:`flyte_bio.modules.chelae` — adapter trimming of FASTQ reads.
- :mod:`flyte_bio.modules.chopper` — long-read quality and length filtering.
- :mod:`flyte_bio.modules.chromap` — fast alignment and preprocessing of chromatin profiles.
- :mod:`flyte_bio.modules.chromsizes` — chromosome sizes and FASTA index.
- :mod:`flyte_bio.modules.clair3` — deep-learning germline small-variant calling from long reads.
- :mod:`flyte_bio.modules.clipkit` — multiple-alignment trimming.
- :mod:`flyte_bio.modules.clustalo` — multiple sequence alignment and guide trees.
- :mod:`flyte_bio.modules.clusty` — Clusty graph clustering.
- :mod:`flyte_bio.modules.cmaple` — maximum-likelihood phylogenies of closely related sequences.
- :mod:`flyte_bio.modules.cnaqc` — CNAqc copy-number QC.
- :mod:`flyte_bio.modules.cnvkit` — CNVkit copy-number calling.
- :mod:`flyte_bio.modules.collectfeaturecounts` — per-sample featureCounts tables merged into one.
- :mod:`flyte_bio.modules.collectstats` — per-sample read statistics joined into one table.
- :mod:`flyte_bio.modules.cooler` — Hi-C contact matrices.
- :mod:`flyte_bio.modules.coreograph` — TMA core detection and cropping.
- :mod:`flyte_bio.modules.coverm` — CoverM read coverage.
- :mod:`flyte_bio.modules.crisprcleanr` — correction of gene-independent CRISPR screen effects.
- :mod:`flyte_bio.modules.csvtk` — table joins and concatenation.
- :mod:`flyte_bio.modules.ctatsplicing` — CTAT-Splicing cancer introns.
- :mod:`flyte_bio.modules.ctree` — clone trees.
- :mod:`flyte_bio.modules.custom` — differential-table and GSEA format helpers.
- :mod:`flyte_bio.modules.cutadapt` — adapter, primer and quality trimming.
- :mod:`flyte_bio.modules.cutesv` — cuteSV SV calling.
- :mod:`flyte_bio.modules.damageprofiler` — ancient-DNA damage patterns.
- :mod:`flyte_bio.modules.datavzrd` — datavzrd HTML reports.
- :mod:`flyte_bio.modules.decoupler` — decoupler activity inference.
- :mod:`flyte_bio.modules.deeptools` — coverage matrices, heatmaps and fingerprint plots.
- :mod:`flyte_bio.modules.deepvariant` — germline small-variant calling.
- :mod:`flyte_bio.modules.demuxem` — demuxEM cell-hashing demultiplexing.
- :mod:`flyte_bio.modules.deseq2` — DESeq2 differential expression.
- :mod:`flyte_bio.modules.dexseq` — differential exon and transcript usage.
- :mod:`flyte_bio.modules.diamond` — protein sequence clustering.
- :mod:`flyte_bio.modules.dragmap` — DRAGMAP short-read alignment (hash table, align).
- :mod:`flyte_bio.modules.dragonflye` — dragonflye bacterial assembly.
- :mod:`flyte_bio.modules.duckdb` — CSV/TSV to Parquet conversion.
- :mod:`flyte_bio.modules.easypqp` — EasyPQP spectral libraries.
- :mod:`flyte_bio.modules.eggnogmapper` — functional annotation of proteins with eggNOG.
- :mod:`flyte_bio.modules.ensemblvep` — variant annotation with Ensembl VEP (cache download, annotate).
- :mod:`flyte_bio.modules.epang` — phylogenetic placement on a reference tree.
- :mod:`flyte_bio.modules.evigene` — transcript set reduction (tr2aacds).
- :mod:`flyte_bio.modules.falco` — falco read QC.
- :mod:`flyte_bio.modules.falint` — fa-lint FASTA validation.
- :mod:`flyte_bio.modules.famsa` — fast multiple sequence alignment and guide trees.
- :mod:`flyte_bio.modules.fastavalidator` — FASTA validation.
- :mod:`flyte_bio.modules.fastplong` — fastplong long-read QC.
- :mod:`flyte_bio.modules.fastqdl` — ENA/SRA FASTQ downloads.
- :mod:`flyte_bio.modules.fastqe` — FASTQ quality summaries as emoji.
- :mod:`flyte_bio.modules.fastqscreen` — screening reads against several genomes.
- :mod:`flyte_bio.modules.fastqutils` — FASTQ validation.
- :mod:`flyte_bio.modules.fcs` — NCBI FCS-adaptor.
- :mod:`flyte_bio.modules.fcsgx` — NCBI FCS-GX contamination screening.
- :mod:`flyte_bio.modules.fibertoolsrs` — m6A and nucleosome calling on Fiber-seq reads.
- :mod:`flyte_bio.modules.filtlong` — long-read quality filtering.
- :mod:`flyte_bio.modules.find` — find and unpigz decompression.
- :mod:`flyte_bio.modules.flye` — long-read de novo assembly.
- :mod:`flyte_bio.modules.fmhfunprofiler` — FracMinHash KEGG ortholog profiling.
- :mod:`flyte_bio.modules.foldmason` — FoldMason structure alignment.
- :mod:`flyte_bio.modules.foldseek` — protein structure database building and search.
- :mod:`flyte_bio.modules.fq` — FASTQ subsampling and linting.
- :mod:`flyte_bio.modules.fqtk` — FASTQ barcode demultiplexing.
- :mod:`flyte_bio.modules.freyja` — Freyja lineage demixing.
- :mod:`flyte_bio.modules.fusioncatcher` — FusionCatcher fusion detection.
- :mod:`flyte_bio.modules.fusioninspector` — FusionInspector fusion validation.
- :mod:`flyte_bio.modules.fusionreport` — fusion-report fusion summaries.
- :mod:`flyte_bio.modules.ganon` — ganon database building.
- :mod:`flyte_bio.modules.gappa` — analysis of phylogenetic placements.
- :mod:`flyte_bio.modules.gapseq` — gapseq metabolic models.
- :mod:`flyte_bio.modules.gatk4` — GATK tools (preprocessing, germline and somatic calling, filtering).
- :mod:`flyte_bio.modules.gawk` — run an awk program over files.
- :mod:`flyte_bio.modules.genmap` — FM index building for mappability.
- :mod:`flyte_bio.modules.genomescope2` — GenomeScope2 genome profiling.
- :mod:`flyte_bio.modules.geoquery` — GEOquery downloads.
- :mod:`flyte_bio.modules.gfaffix` — collapsing of shared affixes in variation graphs.
- :mod:`flyte_bio.modules.gfatools` — gfatools GFA conversion.
- :mod:`flyte_bio.modules.gffread` — GFF/GTF conversion + transcript FASTA.
- :mod:`flyte_bio.modules.glimpse` — GLIMPSE imputation.
- :mod:`flyte_bio.modules.glimpse2` — GLIMPSE2 imputation.
- :mod:`flyte_bio.modules.glnexus` — GLnexus joint genotyping.
- :mod:`flyte_bio.modules.gmmdemux` — GMM-Demux cell-hashing demultiplexing.
- :mod:`flyte_bio.modules.gnusort` — GNU sort of text files.
- :mod:`flyte_bio.modules.gprofiler2` — g:Profiler enrichment.
- :mod:`flyte_bio.modules.graphtyper` — population-scale small-variant genotyping.
- :mod:`flyte_bio.modules.gsea` — Broad GSEA.
- :mod:`flyte_bio.modules.gstama` — TAMA transcript collapse, merge and poly(A) cleanup.
- :mod:`flyte_bio.modules.gtf2bed` — derive a BED12 gene model from a GTF.
- :mod:`flyte_bio.modules.gtffilter` — restrict a GTF to a genome's sequences.
- :mod:`flyte_bio.modules.gtfgenefilter` — keep GTF lines on a genome's sequences (scRNA-seq).
- :mod:`flyte_bio.modules.gunzip` — single-file gzip decompression.
- :mod:`flyte_bio.modules.gzrt` — recover data from corrupt gzipped FASTQs.
- :mod:`flyte_bio.modules.happy` — hap.py benchmarking.
- :mod:`flyte_bio.modules.hasheddrops` — DropletUtils hashedDrops demultiplexing.
- :mod:`flyte_bio.modules.hhsuite` — multiple-alignment format conversion.
- :mod:`flyte_bio.modules.hifiasm` — hifiasm HiFi assembly.
- :mod:`flyte_bio.modules.hificnv` — HiFiCNV (PacBio licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.hiphase` — HiPhase phasing (PacBio licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.hisat2` — splice-aware short-read alignment (splice sites, build, align).
- :mod:`flyte_bio.modules.hmmer` — profile HMM search (hmmsearch).
- :mod:`flyte_bio.modules.homer` — tag directories, peak calling and peak annotation.
- :mod:`flyte_bio.modules.htodemux` — Seurat HTODemux demultiplexing.
- :mod:`flyte_bio.modules.htslib` — bgzip + tabix.
- :mod:`flyte_bio.modules.humann3` — HUMAnN 3 functional profiling.
- :mod:`flyte_bio.modules.humid` — HUMID UMI deduplication.
- :mod:`flyte_bio.modules.ilastik` — pixel classification and multicut segmentation.
- :mod:`flyte_bio.modules.interproscan` — protein domain and family annotation.
- :mod:`flyte_bio.modules.intervals` — scatter-gather genomic intervals.
- :mod:`flyte_bio.modules.iqtree` — maximum-likelihood phylogenies.
- :mod:`flyte_bio.modules.isoseq` — Iso-Seq refine to FLNC reads.
- :mod:`flyte_bio.modules.ivar` — iVar consensus and variants.
- :mod:`flyte_bio.modules.jasminesv` — Jasmine SV merging.
- :mod:`flyte_bio.modules.jellyfish` — Jellyfish k-mer counting.
- :mod:`flyte_bio.modules.kaiju` — Kaiju index building.
- :mod:`flyte_bio.modules.kalign` — Kalign alignment.
- :mod:`flyte_bio.modules.kallisto` — kallisto index build.
- :mod:`flyte_bio.modules.kallistobustools` — kallisto|bustools single-cell counts.
- :mod:`flyte_bio.modules.khmer` — k-mer counting, including unique-k-mer genome size.
- :mod:`flyte_bio.modules.kmcp` — KMCP database building.
- :mod:`flyte_bio.modules.knotannotsv` — knotAnnotSV reports.
- :mod:`flyte_bio.modules.kofamscan` — KEGG Orthology assignment.
- :mod:`flyte_bio.modules.kraken2` — taxonomic classification of reads.
- :mod:`flyte_bio.modules.krakentools` — KrakenTools conversions.
- :mod:`flyte_bio.modules.krakenuniq` — KrakenUniq database building.
- :mod:`flyte_bio.modules.krona` — interactive taxonomy charts.
- :mod:`flyte_bio.modules.last` — LAST genome alignment, training, splitting, conversion and dot plots.
- :mod:`flyte_bio.modules.leafcutter` — intron clustering for differential splicing.
- :mod:`flyte_bio.modules.learnmsa` — learnMSA alignment.
- :mod:`flyte_bio.modules.liftoff` — Liftoff annotation lift-over.
- :mod:`flyte_bio.modules.lima` — PacBio primer and barcode removal.
- :mod:`flyte_bio.modules.limma` — limma differential abundance.
- :mod:`flyte_bio.modules.links` — LINKS scaffolding.
- :mod:`flyte_bio.modules.longphase` — LongPhase phasing.
- :mod:`flyte_bio.modules.longstitch` — LongStitch correction and scaffolding.
- :mod:`flyte_bio.modules.macs2` — ChIP/ATAC-seq peak calling.
- :mod:`flyte_bio.modules.macs3` — model-based peak calling for ChIP-seq and ATAC-seq.
- :mod:`flyte_bio.modules.macsyfinder` — MacSyFinder system detection.
- :mod:`flyte_bio.modules.mafft` — multiple sequence alignment.
- :mod:`flyte_bio.modules.mageck` — CRISPR screen sgRNA counting and gene ranking.
- :mod:`flyte_bio.modules.magus` — MAGUS alignment.
- :mod:`flyte_bio.modules.malt` — MALT index building.
- :mod:`flyte_bio.modules.manta` — Manta somatic SV calling.
- :mod:`flyte_bio.modules.mash` — Mash MinHash sketches.
- :mod:`flyte_bio.modules.mat2json` — MATLAB .mat to JSON/CSV conversion.
- :mod:`flyte_bio.modules.mcquant` — single-cell image quantification.
- :mod:`flyte_bio.modules.md5sum` — MD5 checksums of files.
- :mod:`flyte_bio.modules.mdust` — mdust low-complexity masking.
- :mod:`flyte_bio.modules.medaka` — consensus polishing of nanopore assemblies.
- :mod:`flyte_bio.modules.megahit` — short-read and metagenome assembly.
- :mod:`flyte_bio.modules.megan` — MEGAN Community Edition reports.
- :mod:`flyte_bio.modules.melon` — melon long-read profiling.
- :mod:`flyte_bio.modules.memote` — memote model QC.
- :mod:`flyte_bio.modules.merqury` — Merqury assembly QV.
- :mod:`flyte_bio.modules.meryl` — meryl k-mer databases.
- :mod:`flyte_bio.modules.metacache` — MetaCache database building.
- :mod:`flyte_bio.modules.metaeuk` — eukaryotic gene prediction in contigs.
- :mod:`flyte_bio.modules.metaphlan` — MetaPhlAn taxonomic profiling.
- :mod:`flyte_bio.modules.methyldackel` — MethylDackel methylation calls.
- :mod:`flyte_bio.modules.mgikit` — MGI FASTQ demultiplexing.
- :mod:`flyte_bio.modules.mgnifam` — protein family generation from clustering chunks.
- :mod:`flyte_bio.modules.mifaser` — mi-faser enzyme annotation.
- :mod:`flyte_bio.modules.mindagap` — tile-gap filling and duplicate spot finding.
- :mod:`flyte_bio.modules.minia` — Minia assembly.
- :mod:`flyte_bio.modules.miniasm` — miniasm long-read assembly.
- :mod:`flyte_bio.modules.minimac4` — Minimac4 imputation.
- :mod:`flyte_bio.modules.minimap2` — long- and short-read alignment (index, align).
- :mod:`flyte_bio.modules.mirdeep2` — read mapping and miRNA discovery.
- :mod:`flyte_bio.modules.mirtop` — miRNA annotation in mirGFF3, with counts and export.
- :mod:`flyte_bio.modules.mirtrace` — quality control for small-RNA sequencing.
- :mod:`flyte_bio.modules.misopy` — MISO annotation indexing.
- :mod:`flyte_bio.modules.mmseqs` — MMseqs2 databases, clustering, taxonomy and tables.
- :mod:`flyte_bio.modules.mobster` — MOBSTER subclonal deconvolution.
- :mod:`flyte_bio.modules.modkit` — ONT modified-base calling, pileup and repair (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.mosdepth` — read-depth summaries.
- :mod:`flyte_bio.modules.motus` — mOTUs profiling.
- :mod:`flyte_bio.modules.msisensorpro` — microsatellite site scanning.
- :mod:`flyte_bio.modules.mtmalign` — mTM-align structure alignment (no published licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.multiqc` — aggregate QC report.
- :mod:`flyte_bio.modules.multiqccustombiotype` — biotype counts for MultiQC.
- :mod:`flyte_bio.modules.multiqcsav` — MultiQC reports with Illumina SAV metrics.
- :mod:`flyte_bio.modules.multiseqdemux` — Seurat MULTIseqDemux demultiplexing.
- :mod:`flyte_bio.modules.muscle5` — MUSCLE 5 alignment.
- :mod:`flyte_bio.modules.nacho` — NanoString QC and normalisation.
- :mod:`flyte_bio.modules.nanocomp` — long-read dataset comparison.
- :mod:`flyte_bio.modules.nanolyse` — lambda phage read removal.
- :mod:`flyte_bio.modules.nanoplot` — long-read QC plots and statistics.
- :mod:`flyte_bio.modules.nanoq` — nanoq read filtering.
- :mod:`flyte_bio.modules.ncbigenomedownload` — NCBI genome downloads.
- :mod:`flyte_bio.modules.nextclade` — Nextclade clade assignment.
- :mod:`flyte_bio.modules.nonpareil` — Nonpareil coverage estimates.
- :mod:`flyte_bio.modules.numorph` — light-sheet preprocessing, resampling and nuclei counting.
- :mod:`flyte_bio.modules.odgi` — building, sorting, layout, drawing and stats of variation graphs.
- :mod:`flyte_bio.modules.openms` — OpenMS proteomics tools.
- :mod:`flyte_bio.modules.openmsthirdparty` — Comet search via OpenMS.
- :mod:`flyte_bio.modules.opt` — Off-target probe tracker.
- :mod:`flyte_bio.modules.optitype` — HLA class I genotyping.
- :mod:`flyte_bio.modules.orthofinder` — OrthoFinder orthology.
- :mod:`flyte_bio.modules.pairix` — indexing of contact pairs files.
- :mod:`flyte_bio.modules.pairtools` — Hi-C contact pair parsing, sorting, deduplication and filtering.
- :mod:`flyte_bio.modules.pangolin` — Pangolin lineage assignment.
- :mod:`flyte_bio.modules.parabricks` — NVIDIA Parabricks on GPUs (NVIDIA licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.parsesdrf` — SDRF sample sheet conversion.
- :mod:`flyte_bio.modules.pbccs` — PacBio HiFi consensus reads.
- :mod:`flyte_bio.modules.pbcpgtools` — pb-CpG-tools methylation (PacBio licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.pbjasmine` — 5mC calling on PacBio HiFi reads.
- :mod:`flyte_bio.modules.pbmm2` — minimap2 alignment for PacBio reads.
- :mod:`flyte_bio.modules.pbsv` — pbsv SV calling (PacBio licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.pbtk` — pbtk BAM merging.
- :mod:`flyte_bio.modules.pcgr` — PCGR reference bundle download.
- :mod:`flyte_bio.modules.pear` — paired-end read merging (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.phantompeakqualtools` — ChIP-seq strand cross-correlation QC.
- :mod:`flyte_bio.modules.picard` — Picard MarkDuplicates.
- :mod:`flyte_bio.modules.pigz` — parallel gzip compression and decompression.
- :mod:`flyte_bio.modules.pilon` — Pilon polishing.
- :mod:`flyte_bio.modules.pints` — peak calling for nascent transcription data.
- :mod:`flyte_bio.modules.pirate` — prokaryote pangenome gene families.
- :mod:`flyte_bio.modules.plasmidid` — plasmidID.
- :mod:`flyte_bio.modules.popscle` — popscle pileup, demuxlet and freemuxlet.
- :mod:`flyte_bio.modules.porechop` — Nanopore adapter trimming.
- :mod:`flyte_bio.modules.preseq` — library complexity extrapolation.
- :mod:`flyte_bio.modules.pridepy` — PRIDE archive downloads.
- :mod:`flyte_bio.modules.prinseqplusplus` — PRINSEQ++ read filtering.
- :mod:`flyte_bio.modules.prodigal` — prokaryotic gene prediction.
- :mod:`flyte_bio.modules.prokka` — prokaryotic and viral genome annotation.
- :mod:`flyte_bio.modules.propr` — propr proportionality.
- :mod:`flyte_bio.modules.proteus` — proteus protein groups.
- :mod:`flyte_bio.modules.pyclonevi` — PyClone-VI clustering.
- :mod:`flyte_bio.modules.pycoqc` — pycoQC nanopore QC.
- :mod:`flyte_bio.modules.pydamage` — per-contig ancient-DNA damage estimates.
- :mod:`flyte_bio.modules.qcatch` — cell-level QC of a simpleaf quantification.
- :mod:`flyte_bio.modules.qualimap` — RNA-seq alignment QC.
- :mod:`flyte_bio.modules.quarto` — render parameterised notebooks.
- :mod:`flyte_bio.modules.quast` — assembly quality assessment.
- :mod:`flyte_bio.modules.quilt` — QUILT imputation.
- :mod:`flyte_bio.modules.racon` — consensus polishing of long-read assemblies.
- :mod:`flyte_bio.modules.ragtag` — RagTag scaffolding and patching.
- :mod:`flyte_bio.modules.rastair` — rastair TAPS methylation calls.
- :mod:`flyte_bio.modules.rasusa` — read subsampling to a depth.
- :mod:`flyte_bio.modules.raven` — Raven long-read assembly.
- :mod:`flyte_bio.modules.rclone` — copy and verify data between storage remotes.
- :mod:`flyte_bio.modules.regtools` — splice junction extraction.
- :mod:`flyte_bio.modules.repeatmasker` — RepeatMasker.
- :mod:`flyte_bio.modules.repeatmodeler` — RepeatModeler.
- :mod:`flyte_bio.modules.resolvetaxonomy` — taxonomy tables from taxonomy files or headers.
- :mod:`flyte_bio.modules.rgi` — RGI antimicrobial resistance genes (CARD data licensed; in ``licensed_env``).
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
- :mod:`flyte_bio.modules.samblaster` — samblaster duplicate marking.
- :mod:`flyte_bio.modules.samshee` — Illumina sample sheet validation.
- :mod:`flyte_bio.modules.scanpy` — AnnData conversions (simpleaf to h5ad, concatenation).
- :mod:`flyte_bio.modules.simpleaf` — single-cell quantification with alevin-fry (index, quant).
- :mod:`flyte_bio.modules.samtools` — alignment sort/index/stats/faidx.
- :mod:`flyte_bio.modules.sawfish` — sawfish SV calling (PacBio licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.scimap` — scimap cell clustering.
- :mod:`flyte_bio.modules.seacr` — CUT&RUN/CUT&Tag peak calling.
- :mod:`flyte_bio.modules.sentieon` — Sentieon (untested; licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.seq2hla` — HLA typing and expression from RNA-seq.
- :mod:`flyte_bio.modules.seqcluster` — collapsing of small-RNA reads.
- :mod:`flyte_bio.modules.seqfu` — sequence-length statistics.
- :mod:`flyte_bio.modules.seqkit` — FASTA/FASTQ statistics and pattern replacement.
- :mod:`flyte_bio.modules.seqtk` — FASTQ subsampling.
- :mod:`flyte_bio.modules.sequali` — sequencing-read QC reports.
- :mod:`flyte_bio.modules.seqwish` — variation graph induction from alignments.
- :mod:`flyte_bio.modules.sgdemux` — Singular G4 demultiplexing (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.shapeit5` — SHAPEIT5 phasing.
- :mod:`flyte_bio.modules.shinyngs` — shinyngs reports and apps.
- :mod:`flyte_bio.modules.sigprofiler` — SigProfiler mutational signatures.
- :mod:`flyte_bio.modules.smoothxg` — variation graph smoothing.
- :mod:`flyte_bio.modules.sniffles` — Sniffles2 SV calling.
- :mod:`flyte_bio.modules.snpeff` — variant annotation with snpEff (database download, annotate).
- :mod:`flyte_bio.modules.sortmerna` — rRNA read filtering.
- :mod:`flyte_bio.modules.souporcell` — souporcell genotype-free demultiplexing.
- :mod:`flyte_bio.modules.sourcepredict` — metagenome source prediction.
- :mod:`flyte_bio.modules.sourmash` — MinHash sketches.
- :mod:`flyte_bio.modules.spades` — de novo genome, RNA, viral and metagenome assembly.
- :mod:`flyte_bio.modules.sratools` — SRA downloads and FASTQ dumps.
- :mod:`flyte_bio.modules.spaceranger` — 10x Space Ranger (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.sparsesignatures` — SparseSignatures mutational signatures.
- :mod:`flyte_bio.modules.star` — spliced RNA-seq aligner (index, align).
- :mod:`flyte_bio.modules.stardist` — nucleus segmentation.
- :mod:`flyte_bio.modules.starfusion` — STAR-Fusion fusion detection.
- :mod:`flyte_bio.modules.stitch` — STITCH imputation.
- :mod:`flyte_bio.modules.stranger` — Stranger repeat-expansion annotation.
- :mod:`flyte_bio.modules.strelka` — Strelka somatic small-variant calling.
- :mod:`flyte_bio.modules.stringtie` — transcript assembly / quantification.
- :mod:`flyte_bio.modules.subread` — featureCounts.
- :mod:`flyte_bio.modules.summarizedexperiment` — bundle matrices into an RDS.
- :mod:`flyte_bio.modules.suppa` — alternative splicing events and differential splicing.
- :mod:`flyte_bio.modules.survivor` — SURVIVOR SV utilities.
- :mod:`flyte_bio.modules.svanalyzer` — SVanalyzer benchmarking.
- :mod:`flyte_bio.modules.svim` — SVIM SV calling.
- :mod:`flyte_bio.modules.svtk` — svtk SV standardisation.
- :mod:`flyte_bio.modules.svync` — svync SV harmonisation.
- :mod:`flyte_bio.modules.sylph` — metagenomic profiling by k-mer containment ANI.
- :mod:`flyte_bio.modules.sylphtax` — taxonomic profiles from sylph results.
- :mod:`flyte_bio.modules.taxonkit` — NCBI taxonomy lineages.
- :mod:`flyte_bio.modules.taxpasta` — taxpasta profile standardisation.
- :mod:`flyte_bio.modules.tcoffee` — T-Coffee alignment and scoring.
- :mod:`flyte_bio.modules.thermorawfileparser` — Thermo RAW conversion.
- :mod:`flyte_bio.modules.tiara` — Tiara contig classification.
- :mod:`flyte_bio.modules.tidk` — tidk telomere detection.
- :mod:`flyte_bio.modules.tinc` — TINC tumour-in-normal contamination.
- :mod:`flyte_bio.modules.toulligqc` — Nanopore run QC reports.
- :mod:`flyte_bio.modules.traitar` — Traitar trait prediction.
- :mod:`flyte_bio.modules.transdecoder` — coding regions in transcripts.
- :mod:`flyte_bio.modules.trgt` — TRGT tandem repeats (PacBio licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.trimgalore` — adapter/quality trimming.
- :mod:`flyte_bio.modules.trimmomatic` — Illumina adapter and quality trimming.
- :mod:`flyte_bio.modules.trinity` — de novo transcriptome assembly.
- :mod:`flyte_bio.modules.trnascanse` — tRNAscan-SE tRNA detection.
- :mod:`flyte_bio.modules.trust4` — TCR/BCR repertoire reconstruction.
- :mod:`flyte_bio.modules.truvari` — Truvari SV benchmarking.
- :mod:`flyte_bio.modules.tx2gene` — transcript → gene table from a GTF.
- :mod:`flyte_bio.modules.tximport` — count/TPM matrices from quantifications.
- :mod:`flyte_bio.modules.ucsc` — bedClip / bedGraphToBigWig.
- :mod:`flyte_bio.modules.ultra` — annotation-guided long-read splice alignment.
- :mod:`flyte_bio.modules.umicollapse` — UMI deduplication.
- :mod:`flyte_bio.modules.umitools` — UMI extraction and deduplication.
- :mod:`flyte_bio.modules.unicycler` — Unicycler bacterial assembly.
- :mod:`flyte_bio.modules.untar` — tar archive extraction.
- :mod:`flyte_bio.modules.untarfiles` — tar extraction with a file listing.
- :mod:`flyte_bio.modules.unzip` — zip archive extraction.
- :mod:`flyte_bio.modules.unzipfiles` — archive extraction with a file list.
- :mod:`flyte_bio.modules.upp` — UPP alignment.
- :mod:`flyte_bio.modules.variancepartition` — variancePartition dream.
- :mod:`flyte_bio.modules.variantextractor` — variant-extractor normalisation.
- :mod:`flyte_bio.modules.vcf2maf` — VCF to MAF conversion.
- :mod:`flyte_bio.modules.vcflib` — VCF filtering by field expressions.
- :mod:`flyte_bio.modules.vcftools` — VCF summaries (Ts/Tv, filters).
- :mod:`flyte_bio.modules.viber` — VIBER read-count clustering.
- :mod:`flyte_bio.modules.vireo` — vireo donor assignment.
- :mod:`flyte_bio.modules.vrhyme` — vRhyme viral binning.
- :mod:`flyte_bio.modules.vsearch` — sequence clustering.
- :mod:`flyte_bio.modules.vt` — vt variant decomposition and normalisation.
- :mod:`flyte_bio.modules.wfmash` — all-to-all genome alignment for pangenome graphs.
- :mod:`flyte_bio.modules.wget` — file download with GNU Wget.
- :mod:`flyte_bio.modules.wipertools` — repair malformed FASTQs in parallel chunks.
- :mod:`flyte_bio.modules.wittyer` — Illumina witty.er (Illumina licence; in ``licensed_env``).
- :mod:`flyte_bio.modules.xeniumranger` — 10x Xenium Ranger (licensed; in ``licensed_env``).
- :mod:`flyte_bio.modules.xz` — xz decompression.
- :mod:`flyte_bio.modules.yahs` — YaHS Hi-C scaffolding.
- :mod:`flyte_bio.modules.yara` — read mapping (index, mapper).
- :mod:`flyte_bio.modules.zip` — 7-Zip archives.

The module-level :data:`env` here is an aggregate
:class:`flyte.TaskEnvironment` depending on every submodule's env.
Pipelines depend on it once to gain access to every wrapped tool::

    from flyte_bio.modules import env as modules_env
"""

import flyte

from .aardvark import env as aardvark_env
from .abacas import env as abacas_env
from .adapterremoval import env as adapterremoval_env
from .affy import env as affy_env
from .agat import env as agat_env
from .amulety import env as amulety_env
from .anndata import env as anndata_env
from .anndatar import env as anndatar_env
from .annotsv import env as annotsv_env
from .aria2 import env as aria2_env
from .arriba import env as arriba_env
from .art import env as art_env
from .artic import env as artic_env
from .ashlar import env as ashlar_env
from .assemblyscan import env as assemblyscan_env
from .ataqv import env as ataqv_env
from .atlasgeneannotationmanipulation import env as atlasgeneannotationmanipulation_env
from .autocycler import env as autocycler_env
from .backsub import env as backsub_env
from .bakta import env as bakta_env
from .bamtools import env as bamtools_env
from .bandage import env as bandage_env
from .barrnap import env as barrnap_env
from .bases2fastq import env as bases2fastq_env
from .basicpy import env as basicpy_env
from .bbmap import env as bbmap_env
from .bcftools import env as bcftools_env
from .bcl2fastq import env as bcl2fastq_env
from .bclconvert import env as bclconvert_env
from .beagle5 import env as beagle5_env
from .bedops import env as bedops_env
from .bedtools import env as bedtools_env
from .bedtools import genomecov_env
from .bff import env as bff_env
from .bftools import env as bftools_env
from .bioawk import env as bioawk_env
from .bismark import env as bismark_env
from .blast import env as blast_env
from .blat import env as blat_env
from .bowtie import env as bowtie_env
from .bowtie2 import env as bowtie2_env
from .bracken import env as bracken_env
from .busco import env as busco_env
from .bwa import env as bwa_env
from .bwamem2 import env as bwamem2_env
from .bwameth import env as bwameth_env
from .canu import env as canu_env
from .carveme import env as carveme_env
from .cat import env as cat_env
from .catadditionalfasta import env as catadditionalfasta_env
from .catpack import env as catpack_env
from .cdhit import env as cdhit_env
from .cellbender import env as cellbender_env
from .cellpose import env as cellpose_env
from .cellranger import env as cellranger_env
from .cellrangerarc import env as cellrangerarc_env
from .cellsnp import env as cellsnp_env
from .centrifuge import env as centrifuge_env
from .centrifuger import env as centrifuger_env
from .checkm2 import env as checkm2_env
from .checkqc import env as checkqc_env
from .checkv import env as checkv_env
from .chelae import env as chelae_env
from .chopper import env as chopper_env
from .chromap import env as chromap_env
from .chromsizes import env as chromsizes_env
from .clair3 import env as clair3_env
from .clipkit import env as clipkit_env
from .clustalo import env as clustalo_env
from .clusty import env as clusty_env
from .cmaple import env as cmaple_env
from .cnaqc import env as cnaqc_env
from .cnvkit import env as cnvkit_env
from .collectfeaturecounts import env as collectfeaturecounts_env
from .collectstats import env as collectstats_env
from .cooler import env as cooler_env
from .coreograph import env as coreograph_env
from .coverm import env as coverm_env
from .crisprcleanr import env as crisprcleanr_env
from .csvtk import env as csvtk_env
from .ctatsplicing import env as ctatsplicing_env
from .ctree import env as ctree_env
from .custom import env as custom_env
from .cutadapt import env as cutadapt_env
from .cutesv import env as cutesv_env
from .damageprofiler import env as damageprofiler_env
from .datavzrd import env as datavzrd_env
from .decoupler import env as decoupler_env
from .deepcell import env as deepcell_env
from .deeptools import env as deeptools_env
from .deepvariant import env as deepvariant_env
from .demuxem import env as demuxem_env
from .deseq2 import env as deseq2_env
from .deseq2_qc import env as deseq2_qc_env
from .dexseq import env as dexseq_env
from .diamond import env as diamond_env
from .dragmap import env as dragmap_env
from .dragonflye import env as dragonflye_env
from .duckdb import env as duckdb_env
from .dupradar import env as dupradar_env
from .easypqp import env as easypqp_env
from .eggnogmapper import env as eggnogmapper_env
from .ensemblvep import env as ensemblvep_env
from .epang import env as epang_env
from .evigene import env as evigene_env
from .falco import env as falco_env
from .falint import env as falint_env
from .famsa import env as famsa_env
from .fastavalidator import env as fastavalidator_env
from .fastp import env as fastp_env
from .fastplong import env as fastplong_env
from .fastqc import env as fastqc_env
from .fastqdl import env as fastqdl_env
from .fastqe import env as fastqe_env
from .fastqscreen import env as fastqscreen_env
from .fastqutils import env as fastqutils_env
from .fcs import env as fcs_env
from .fcsgx import env as fcsgx_env
from .fibertoolsrs import env as fibertoolsrs_env
from .filtlong import env as filtlong_env
from .find import env as find_env
from .flye import env as flye_env
from .fmhfunprofiler import env as fmhfunprofiler_env
from .foldmason import env as foldmason_env
from .foldseek import env as foldseek_env
from .fq import env as fq_env
from .fqtk import env as fqtk_env
from .freyja import env as freyja_env
from .fusioncatcher import env as fusioncatcher_env
from .fusioninspector import env as fusioninspector_env
from .fusionreport import env as fusionreport_env
from .ganon import env as ganon_env
from .gappa import env as gappa_env
from .gapseq import env as gapseq_env
from .gatk4 import env as gatk4_env
from .gawk import env as gawk_env
from .genmap import env as genmap_env
from .genomescope2 import env as genomescope2_env
from .geoquery import env as geoquery_env
from .gfaffix import env as gfaffix_env
from .gfatools import env as gfatools_env
from .gffread import env as gffread_env
from .glimpse import env as glimpse_env
from .glimpse2 import env as glimpse2_env
from .glnexus import env as glnexus_env
from .gmmdemux import env as gmmdemux_env
from .gnusort import env as gnusort_env
from .gprofiler2 import env as gprofiler2_env
from .graphtyper import env as graphtyper_env
from .gsea import env as gsea_env
from .gstama import env as gstama_env
from .gtf2bed import env as gtf2bed_env
from .gtffilter import env as gtffilter_env
from .gtfgenefilter import env as gtfgenefilter_env
from .gunzip import env as gunzip_env
from .gzrt import env as gzrt_env
from .happy import env as happy_env
from .hasheddrops import env as hasheddrops_env
from .hhsuite import env as hhsuite_env
from .hifiasm import env as hifiasm_env
from .hificnv import env as hificnv_env
from .hiphase import env as hiphase_env
from .hisat2 import env as hisat2_env
from .hmmer import env as hmmer_env
from .homer import env as homer_env
from .htodemux import env as htodemux_env
from .htslib import env as htslib_env
from .humann3 import env as humann3_env
from .humid import env as humid_env
from .ilastik import env as ilastik_env
from .interproscan import env as interproscan_env
from .intervals import env as intervals_env
from .iqtree import env as iqtree_env
from .isoseq import env as isoseq_env
from .ivar import env as ivar_env
from .jasminesv import env as jasminesv_env
from .jellyfish import env as jellyfish_env
from .kaiju import env as kaiju_env
from .kalign import env as kalign_env
from .kallisto import env as kallisto_env
from .kallistobustools import env as kallistobustools_env
from .khmer import env as khmer_env
from .kmcp import env as kmcp_env
from .knotannotsv import env as knotannotsv_env
from .kofamscan import env as kofamscan_env
from .kraken2 import env as kraken2_env
from .krakentools import env as krakentools_env
from .krakenuniq import env as krakenuniq_env
from .krona import env as krona_env
from .last import env as last_env
from .leafcutter import env as leafcutter_env
from .learnmsa import env as learnmsa_env
from .liftoff import env as liftoff_env
from .lima import env as lima_env
from .limma import env as limma_env
from .links import env as links_env
from .longphase import env as longphase_env
from .longstitch import env as longstitch_env
from .macs2 import env as macs2_env
from .macs3 import env as macs3_env
from .macsyfinder import env as macsyfinder_env
from .mafft import env as mafft_env
from .mageck import env as mageck_env
from .magus import env as magus_env
from .malt import env as malt_env
from .manta import env as manta_env
from .mash import env as mash_env
from .mat2json import env as mat2json_env
from .mcquant import env as mcquant_env
from .md5sum import env as md5sum_env
from .mdust import env as mdust_env
from .medaka import env as medaka_env
from .megahit import env as megahit_env
from .megan import env as megan_env
from .melon import env as melon_env
from .memote import env as memote_env
from .merqury import env as merqury_env
from .meryl import env as meryl_env
from .metacache import env as metacache_env
from .metaeuk import env as metaeuk_env
from .metaphlan import env as metaphlan_env
from .methyldackel import env as methyldackel_env
from .mgikit import env as mgikit_env
from .mgnifam import env as mgnifam_env
from .mifaser import env as mifaser_env
from .mindagap import env as mindagap_env
from .minia import env as minia_env
from .miniasm import env as miniasm_env
from .minimac4 import env as minimac4_env
from .minimap2 import env as minimap2_env
from .mirdeep2 import env as mirdeep2_env
from .mirtop import env as mirtop_env
from .mirtrace import env as mirtrace_env
from .misopy import env as misopy_env
from .mmseqs import env as mmseqs_env
from .mobster import env as mobster_env
from .modkit import env as modkit_env
from .mosdepth import env as mosdepth_env
from .motus import env as motus_env
from .msisensorpro import env as msisensorpro_env
from .mtmalign import env as mtmalign_env
from .multiqc import env as multiqc_env
from .multiqccustombiotype import env as multiqccustombiotype_env
from .multiqcsav import env as multiqcsav_env
from .multiseqdemux import env as multiseqdemux_env
from .muscle5 import env as muscle5_env
from .nacho import env as nacho_env
from .nanocomp import env as nanocomp_env
from .nanolyse import env as nanolyse_env
from .nanoplot import env as nanoplot_env
from .nanoq import env as nanoq_env
from .ncbigenomedownload import env as ncbigenomedownload_env
from .nextclade import env as nextclade_env
from .nonpareil import env as nonpareil_env
from .numorph import env as numorph_env
from .odgi import env as odgi_env
from .openms import env as openms_env
from .openmsthirdparty import env as openmsthirdparty_env
from .opt import env as opt_env
from .optitype import env as optitype_env
from .orthofinder import env as orthofinder_env
from .pairix import env as pairix_env
from .pairtools import env as pairtools_env
from .pangolin import env as pangolin_env
from .parabricks import env as parabricks_env
from .parsesdrf import env as parsesdrf_env
from .pbccs import env as pbccs_env
from .pbcpgtools import env as pbcpgtools_env
from .pbjasmine import env as pbjasmine_env
from .pbmm2 import env as pbmm2_env
from .pbsv import env as pbsv_env
from .pbtk import env as pbtk_env
from .pcgr import env as pcgr_env
from .pear import env as pear_env
from .phantompeakqualtools import env as phantompeakqualtools_env
from .picard import env as picard_env
from .pigz import env as pigz_env
from .pilon import env as pilon_env
from .pints import env as pints_env
from .pirate import env as pirate_env
from .plasmidid import env as plasmidid_env
from .popscle import env as popscle_env
from .porechop import env as porechop_env
from .preseq import env as preseq_env
from .pridepy import env as pridepy_env
from .prinseqplusplus import env as prinseqplusplus_env
from .prodigal import env as prodigal_env
from .prokka import env as prokka_env
from .propr import env as propr_env
from .proteus import env as proteus_env
from .pyclonevi import env as pyclonevi_env
from .pycoqc import env as pycoqc_env
from .pydamage import env as pydamage_env
from .qcatch import env as qcatch_env
from .qualimap import env as qualimap_env
from .quarto import env as quarto_env
from .quast import env as quast_env
from .quilt import env as quilt_env
from .racon import env as racon_env
from .ragtag import env as ragtag_env
from .rastair import env as rastair_env
from .rasusa import env as rasusa_env
from .raven import env as raven_env
from .rclone import env as rclone_env
from .regtools import env as regtools_env
from .repeatmasker import env as repeatmasker_env
from .repeatmodeler import env as repeatmodeler_env
from .resolvetaxonomy import env as resolvetaxonomy_env
from .rgi import env as rgi_env
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
from .samblaster import env as samblaster_env
from .samshee import env as samshee_env
from .samtools import env as samtools_env
from .sawfish import env as sawfish_env
from .scanpy import env as scanpy_env
from .scimap import env as scimap_env
from .seacr import env as seacr_env
from .sentieon import env as sentieon_env
from .seq2hla import env as seq2hla_env
from .seqcluster import env as seqcluster_env
from .seqfu import env as seqfu_env
from .seqkit import env as seqkit_env
from .seqtk import env as seqtk_env
from .sequali import env as sequali_env
from .seqwish import env as seqwish_env
from .sgdemux import env as sgdemux_env
from .shapeit5 import env as shapeit5_env
from .shinyngs import env as shinyngs_env
from .sigprofiler import env as sigprofiler_env
from .simpleaf import env as simpleaf_env
from .smoothxg import env as smoothxg_env
from .sniffles import env as sniffles_env
from .snpeff import env as snpeff_env
from .sortmerna import env as sortmerna_env
from .souporcell import env as souporcell_env
from .sourcepredict import env as sourcepredict_env
from .sourmash import env as sourmash_env
from .spaceranger import env as spaceranger_env
from .spades import env as spades_env
from .sparsesignatures import env as sparsesignatures_env
from .sratools import env as sratools_env
from .star import env as star_env
from .stardist import env as stardist_env
from .starfusion import env as starfusion_env
from .stitch import env as stitch_env
from .stranger import env as stranger_env
from .strelka import env as strelka_env
from .stringtie import env as stringtie_env
from .subread import env as subread_env
from .summarizedexperiment import env as summarizedexperiment_env
from .suppa import env as suppa_env
from .survivor import env as survivor_env
from .svanalyzer import env as svanalyzer_env
from .svim import env as svim_env
from .svtk import env as svtk_env
from .svync import env as svync_env
from .sylph import env as sylph_env
from .sylphtax import env as sylphtax_env
from .taxonkit import env as taxonkit_env
from .taxpasta import env as taxpasta_env
from .tcoffee import env as tcoffee_env
from .thermorawfileparser import env as thermorawfileparser_env
from .tiara import env as tiara_env
from .tidk import env as tidk_env
from .tinc import env as tinc_env
from .toulligqc import env as toulligqc_env
from .traitar import env as traitar_env
from .transdecoder import env as transdecoder_env
from .trgt import env as trgt_env
from .trimgalore import env as trimgalore_env
from .trimmomatic import env as trimmomatic_env
from .trinity import env as trinity_env
from .trnascanse import env as trnascanse_env
from .trust4 import env as trust4_env
from .truvari import env as truvari_env
from .tx2gene import env as tx2gene_env
from .tximport import env as tximport_env
from .ucsc import env as ucsc_env
from .ultra import env as ultra_env
from .umicollapse import env as umicollapse_env
from .umitools import env as umitools_env
from .unicycler import env as unicycler_env
from .untar import env as untar_env
from .untarfiles import env as untarfiles_env
from .unzip import env as unzip_env
from .unzipfiles import env as unzipfiles_env
from .upp import env as upp_env
from .variancepartition import env as variancepartition_env
from .variantextractor import env as variantextractor_env
from .vcf2maf import env as vcf2maf_env
from .vcflib import env as vcflib_env
from .vcftools import env as vcftools_env
from .viber import env as viber_env
from .vireo import env as vireo_env
from .vrhyme import env as vrhyme_env
from .vsearch import env as vsearch_env
from .vt import env as vt_env
from .wfmash import env as wfmash_env
from .wget import env as wget_env
from .wipertools import env as wipertools_env
from .wittyer import env as wittyer_env
from .xeniumranger import env as xeniumranger_env
from .xz import env as xz_env
from .yahs import env as yahs_env
from .yara import env as yara_env
from .zip import env as zip_env

env = flyte.TaskEnvironment(
    name="flyte_bio_modules",
    depends_on=[
        aardvark_env,
        abacas_env,
        adapterremoval_env,
        affy_env,
        agat_env,
        amulety_env,
        anndata_env,
        anndatar_env,
        aria2_env,
        arriba_env,
        art_env,
        artic_env,
        ashlar_env,
        assemblyscan_env,
        ataqv_env,
        atlasgeneannotationmanipulation_env,
        autocycler_env,
        backsub_env,
        bakta_env,
        bamtools_env,
        bandage_env,
        barrnap_env,
        basicpy_env,
        bbmap_env,
        bcftools_env,
        beagle5_env,
        bedops_env,
        bedtools_env,
        bff_env,
        bftools_env,
        bioawk_env,
        bismark_env,
        blast_env,
        bowtie2_env,
        bowtie_env,
        bracken_env,
        busco_env,
        bwa_env,
        bwamem2_env,
        bwameth_env,
        canu_env,
        carveme_env,
        cat_env,
        catadditionalfasta_env,
        catpack_env,
        cdhit_env,
        cellbender_env,
        cellpose_env,
        cellsnp_env,
        centrifuge_env,
        centrifuger_env,
        checkm2_env,
        checkqc_env,
        checkv_env,
        chelae_env,
        chopper_env,
        chromap_env,
        chromsizes_env,
        clair3_env,
        clipkit_env,
        clustalo_env,
        clusty_env,
        cmaple_env,
        cnaqc_env,
        cnvkit_env,
        collectfeaturecounts_env,
        collectstats_env,
        cooler_env,
        coreograph_env,
        coverm_env,
        crisprcleanr_env,
        csvtk_env,
        ctatsplicing_env,
        ctree_env,
        custom_env,
        cutadapt_env,
        cutesv_env,
        damageprofiler_env,
        datavzrd_env,
        decoupler_env,
        deeptools_env,
        deepvariant_env,
        demuxem_env,
        deseq2_env,
        deseq2_qc_env,
        dexseq_env,
        diamond_env,
        dragmap_env,
        dragonflye_env,
        duckdb_env,
        dupradar_env,
        easypqp_env,
        eggnogmapper_env,
        ensemblvep_env,
        epang_env,
        evigene_env,
        falco_env,
        falint_env,
        famsa_env,
        fastavalidator_env,
        fastp_env,
        fastplong_env,
        fastqc_env,
        fastqdl_env,
        fastqe_env,
        fastqscreen_env,
        fastqutils_env,
        fcs_env,
        fcsgx_env,
        fibertoolsrs_env,
        filtlong_env,
        find_env,
        flye_env,
        fmhfunprofiler_env,
        foldmason_env,
        foldseek_env,
        fq_env,
        fqtk_env,
        freyja_env,
        fusioncatcher_env,
        fusioninspector_env,
        fusionreport_env,
        ganon_env,
        gappa_env,
        gapseq_env,
        gatk4_env,
        gawk_env,
        genmap_env,
        genomecov_env,
        genomescope2_env,
        geoquery_env,
        gfaffix_env,
        gfatools_env,
        gffread_env,
        glimpse2_env,
        glimpse_env,
        glnexus_env,
        gmmdemux_env,
        gnusort_env,
        gprofiler2_env,
        graphtyper_env,
        gsea_env,
        gstama_env,
        gtf2bed_env,
        gtffilter_env,
        gtfgenefilter_env,
        gunzip_env,
        gzrt_env,
        happy_env,
        hasheddrops_env,
        hhsuite_env,
        hifiasm_env,
        hisat2_env,
        hmmer_env,
        homer_env,
        htodemux_env,
        htslib_env,
        humann3_env,
        humid_env,
        ilastik_env,
        interproscan_env,
        intervals_env,
        iqtree_env,
        isoseq_env,
        ivar_env,
        jasminesv_env,
        jellyfish_env,
        kaiju_env,
        kalign_env,
        kallisto_env,
        kallistobustools_env,
        khmer_env,
        kmcp_env,
        knotannotsv_env,
        kofamscan_env,
        kraken2_env,
        krakentools_env,
        krakenuniq_env,
        krona_env,
        last_env,
        leafcutter_env,
        learnmsa_env,
        liftoff_env,
        lima_env,
        limma_env,
        links_env,
        longphase_env,
        longstitch_env,
        macs2_env,
        macs3_env,
        macsyfinder_env,
        mafft_env,
        mageck_env,
        magus_env,
        malt_env,
        manta_env,
        mash_env,
        mat2json_env,
        mcquant_env,
        md5sum_env,
        mdust_env,
        medaka_env,
        megahit_env,
        megan_env,
        melon_env,
        memote_env,
        merqury_env,
        meryl_env,
        metacache_env,
        metaeuk_env,
        metaphlan_env,
        methyldackel_env,
        mgikit_env,
        mgnifam_env,
        mifaser_env,
        mindagap_env,
        minia_env,
        miniasm_env,
        minimac4_env,
        minimap2_env,
        mirdeep2_env,
        mirtop_env,
        mirtrace_env,
        misopy_env,
        mmseqs_env,
        mobster_env,
        mosdepth_env,
        motus_env,
        msisensorpro_env,
        multiqc_env,
        multiqccustombiotype_env,
        multiqcsav_env,
        multiseqdemux_env,
        muscle5_env,
        nacho_env,
        nanocomp_env,
        nanolyse_env,
        nanoplot_env,
        nanoq_env,
        ncbigenomedownload_env,
        nextclade_env,
        nonpareil_env,
        numorph_env,
        odgi_env,
        openms_env,
        openmsthirdparty_env,
        opt_env,
        optitype_env,
        orthofinder_env,
        pairix_env,
        pairtools_env,
        pangolin_env,
        parsesdrf_env,
        pbccs_env,
        pbjasmine_env,
        pbmm2_env,
        pbtk_env,
        pcgr_env,
        phantompeakqualtools_env,
        picard_env,
        pigz_env,
        pilon_env,
        pints_env,
        pirate_env,
        plasmidid_env,
        popscle_env,
        porechop_env,
        preseq_env,
        pridepy_env,
        prinseqplusplus_env,
        prodigal_env,
        prokka_env,
        propr_env,
        proteus_env,
        pyclonevi_env,
        pycoqc_env,
        pydamage_env,
        qcatch_env,
        qualimap_env,
        quarto_env,
        quast_env,
        quilt_env,
        racon_env,
        ragtag_env,
        rastair_env,
        rasusa_env,
        raven_env,
        rclone_env,
        regtools_env,
        repeatmasker_env,
        repeatmodeler_env,
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
        samblaster_env,
        samshee_env,
        samtools_env,
        scanpy_env,
        scimap_env,
        seacr_env,
        seq2hla_env,
        seqcluster_env,
        seqfu_env,
        seqkit_env,
        seqtk_env,
        sequali_env,
        seqwish_env,
        shapeit5_env,
        shinyngs_env,
        sigprofiler_env,
        simpleaf_env,
        smoothxg_env,
        sniffles_env,
        snpeff_env,
        sortmerna_env,
        souporcell_env,
        sourcepredict_env,
        sourmash_env,
        spades_env,
        sparsesignatures_env,
        sratools_env,
        star_env,
        stardist_env,
        starfusion_env,
        stitch_env,
        stranger_env,
        strelka_env,
        stringtie_env,
        subread_env,
        summarizedexperiment_env,
        suppa_env,
        survivor_env,
        svanalyzer_env,
        svim_env,
        svtk_env,
        svync_env,
        sylph_env,
        sylphtax_env,
        taxonkit_env,
        taxpasta_env,
        tcoffee_env,
        thermorawfileparser_env,
        tiara_env,
        tidk_env,
        tinc_env,
        toulligqc_env,
        traitar_env,
        transdecoder_env,
        trimgalore_env,
        trimmomatic_env,
        trinity_env,
        trnascanse_env,
        trust4_env,
        truvari_env,
        tx2gene_env,
        tximport_env,
        ucsc_env,
        ultra_env,
        umicollapse_env,
        umitools_env,
        unicycler_env,
        untar_env,
        untarfiles_env,
        unzip_env,
        unzipfiles_env,
        upp_env,
        variancepartition_env,
        variantextractor_env,
        vcf2maf_env,
        vcflib_env,
        vcftools_env,
        viber_env,
        vireo_env,
        vrhyme_env,
        vsearch_env,
        vt_env,
        wfmash_env,
        wget_env,
        wipertools_env,
        xz_env,
        yahs_env,
        yara_env,
        zip_env,
    ],
)

# Tools whose licences restrict use, left out of the default `env`:
# - Cell Ranger and Space Ranger: 10x Genomics' end-user licence (internal use
#   with 10x products or data, no redistribution).
# - DeepCell's Mesmer, PEAR and rMATS: non-commercial use only.
# - BLAT: free for academic, non-profit and personal use; commercial use needs a
#   Kent Informatics licence.
# - modkit: Oxford Nanopore's public licence (research use).
# - bcl2fastq and BCL Convert: Illumina's end-user licence.
# - bases2fastq: Element Biosciences' end-user licence.
# - sgdemux: licensed for the Singular Genomics G4 platform only.
# - HiFiCNV, TRGT, HiPhase, pb-CpG-tools, pbsv and sawfish: the PacBio Software
#   License Agreement (PacBio data only; binary-only; confidentiality).
# - Xenium Ranger and Cell Ranger ARC: 10x Genomics' end-user licence.
# - Parabricks: NVIDIA AI Product Agreement (click-through, no redistribution).
# - witty.er: Illumina licence agreement (not published).
# - Sentieon: commercial; runs only against a licence server (untested here).
# - AnnotSV: its annotation bundle includes OMIM- and COSMIC-derived data that
#   restrict commercial use.
# - mTM-align: published without a licence (all rights reserved).
# - RGI: the CARD database it uses bars commercial organisations without a
#   McMaster licence.
# Depend on `licensed_env` explicitly once you've checked the terms.
licensed_env = flyte.TaskEnvironment(
    name="flyte_bio_licensed_modules",
    depends_on=[
        annotsv_env,
        bases2fastq_env,
        bcl2fastq_env,
        bclconvert_env,
        blat_env,
        cellranger_env,
        cellrangerarc_env,
        deepcell_env,
        hificnv_env,
        hiphase_env,
        modkit_env,
        mtmalign_env,
        parabricks_env,
        pbcpgtools_env,
        pbsv_env,
        pear_env,
        rgi_env,
        rmats_env,
        sawfish_env,
        sentieon_env,
        sgdemux_env,
        spaceranger_env,
        trgt_env,
        wittyer_env,
        xeniumranger_env,
    ],
)

__all__ = [
    "aardvark_env",
    "abacas_env",
    "adapterremoval_env",
    "affy_env",
    "agat_env",
    "amulety_env",
    "anndata_env",
    "anndatar_env",
    "annotsv_env",
    "aria2_env",
    "arriba_env",
    "art_env",
    "artic_env",
    "ashlar_env",
    "assemblyscan_env",
    "ataqv_env",
    "atlasgeneannotationmanipulation_env",
    "autocycler_env",
    "backsub_env",
    "bakta_env",
    "bamtools_env",
    "bandage_env",
    "barrnap_env",
    "bases2fastq_env",
    "basicpy_env",
    "bbmap_env",
    "bcftools_env",
    "bcl2fastq_env",
    "bclconvert_env",
    "beagle5_env",
    "bedops_env",
    "bedtools_env",
    "bff_env",
    "bftools_env",
    "bioawk_env",
    "bismark_env",
    "blast_env",
    "blat_env",
    "bowtie2_env",
    "bowtie_env",
    "bracken_env",
    "busco_env",
    "bwa_env",
    "bwamem2_env",
    "bwameth_env",
    "canu_env",
    "carveme_env",
    "cat_env",
    "catadditionalfasta_env",
    "catpack_env",
    "cdhit_env",
    "cellbender_env",
    "cellpose_env",
    "cellranger_env",
    "cellrangerarc_env",
    "cellsnp_env",
    "centrifuge_env",
    "centrifuger_env",
    "checkm2_env",
    "checkqc_env",
    "checkv_env",
    "chelae_env",
    "chopper_env",
    "chromap_env",
    "chromsizes_env",
    "clair3_env",
    "clipkit_env",
    "clustalo_env",
    "clusty_env",
    "cmaple_env",
    "cnaqc_env",
    "cnvkit_env",
    "collectfeaturecounts_env",
    "collectstats_env",
    "cooler_env",
    "coreograph_env",
    "coverm_env",
    "crisprcleanr_env",
    "csvtk_env",
    "ctatsplicing_env",
    "ctree_env",
    "custom_env",
    "cutadapt_env",
    "cutesv_env",
    "damageprofiler_env",
    "datavzrd_env",
    "decoupler_env",
    "deepcell_env",
    "deeptools_env",
    "deepvariant_env",
    "demuxem_env",
    "deseq2_env",
    "deseq2_qc_env",
    "dexseq_env",
    "diamond_env",
    "dragmap_env",
    "dragonflye_env",
    "duckdb_env",
    "dupradar_env",
    "easypqp_env",
    "eggnogmapper_env",
    "ensemblvep_env",
    "env",
    "epang_env",
    "evigene_env",
    "falco_env",
    "falint_env",
    "famsa_env",
    "fastavalidator_env",
    "fastp_env",
    "fastplong_env",
    "fastqc_env",
    "fastqdl_env",
    "fastqe_env",
    "fastqscreen_env",
    "fastqutils_env",
    "fcs_env",
    "fcsgx_env",
    "fibertoolsrs_env",
    "filtlong_env",
    "find_env",
    "flye_env",
    "fmhfunprofiler_env",
    "foldmason_env",
    "foldseek_env",
    "fq_env",
    "fqtk_env",
    "freyja_env",
    "fusioncatcher_env",
    "fusioninspector_env",
    "fusionreport_env",
    "ganon_env",
    "gappa_env",
    "gapseq_env",
    "gatk4_env",
    "gawk_env",
    "genmap_env",
    "genomecov_env",
    "genomescope2_env",
    "geoquery_env",
    "gfaffix_env",
    "gfatools_env",
    "gffread_env",
    "glimpse2_env",
    "glimpse_env",
    "glnexus_env",
    "gmmdemux_env",
    "gnusort_env",
    "gprofiler2_env",
    "graphtyper_env",
    "gsea_env",
    "gstama_env",
    "gtf2bed_env",
    "gtffilter_env",
    "gtfgenefilter_env",
    "gunzip_env",
    "gzrt_env",
    "happy_env",
    "hasheddrops_env",
    "hhsuite_env",
    "hifiasm_env",
    "hificnv_env",
    "hiphase_env",
    "hisat2_env",
    "hmmer_env",
    "homer_env",
    "htodemux_env",
    "htslib_env",
    "humann3_env",
    "humid_env",
    "ilastik_env",
    "interproscan_env",
    "intervals_env",
    "iqtree_env",
    "isoseq_env",
    "ivar_env",
    "jasminesv_env",
    "jellyfish_env",
    "kaiju_env",
    "kalign_env",
    "kallisto_env",
    "kallistobustools_env",
    "khmer_env",
    "kmcp_env",
    "knotannotsv_env",
    "kofamscan_env",
    "kraken2_env",
    "krakentools_env",
    "krakenuniq_env",
    "krona_env",
    "last_env",
    "leafcutter_env",
    "learnmsa_env",
    "licensed_env",
    "liftoff_env",
    "lima_env",
    "limma_env",
    "links_env",
    "longphase_env",
    "longstitch_env",
    "macs2_env",
    "macs3_env",
    "macsyfinder_env",
    "mafft_env",
    "mageck_env",
    "magus_env",
    "malt_env",
    "manta_env",
    "mash_env",
    "mat2json_env",
    "mcquant_env",
    "md5sum_env",
    "mdust_env",
    "medaka_env",
    "megahit_env",
    "megan_env",
    "melon_env",
    "memote_env",
    "merqury_env",
    "meryl_env",
    "metacache_env",
    "metaeuk_env",
    "metaphlan_env",
    "methyldackel_env",
    "mgikit_env",
    "mgnifam_env",
    "mifaser_env",
    "mindagap_env",
    "minia_env",
    "miniasm_env",
    "minimac4_env",
    "minimap2_env",
    "mirdeep2_env",
    "mirtop_env",
    "mirtrace_env",
    "misopy_env",
    "mmseqs_env",
    "mobster_env",
    "modkit_env",
    "mosdepth_env",
    "motus_env",
    "msisensorpro_env",
    "mtmalign_env",
    "multiqc_env",
    "multiqccustombiotype_env",
    "multiqcsav_env",
    "multiseqdemux_env",
    "muscle5_env",
    "nacho_env",
    "nanocomp_env",
    "nanolyse_env",
    "nanoplot_env",
    "nanoq_env",
    "ncbigenomedownload_env",
    "nextclade_env",
    "nonpareil_env",
    "numorph_env",
    "odgi_env",
    "openms_env",
    "openmsthirdparty_env",
    "opt_env",
    "optitype_env",
    "orthofinder_env",
    "pairix_env",
    "pairtools_env",
    "pangolin_env",
    "parabricks_env",
    "parsesdrf_env",
    "pbccs_env",
    "pbcpgtools_env",
    "pbjasmine_env",
    "pbmm2_env",
    "pbsv_env",
    "pbtk_env",
    "pcgr_env",
    "pear_env",
    "phantompeakqualtools_env",
    "picard_env",
    "pigz_env",
    "pilon_env",
    "pints_env",
    "pirate_env",
    "plasmidid_env",
    "popscle_env",
    "porechop_env",
    "preseq_env",
    "pridepy_env",
    "prinseqplusplus_env",
    "prodigal_env",
    "prokka_env",
    "propr_env",
    "proteus_env",
    "pyclonevi_env",
    "pycoqc_env",
    "pydamage_env",
    "qcatch_env",
    "qualimap_env",
    "quarto_env",
    "quast_env",
    "quilt_env",
    "racon_env",
    "ragtag_env",
    "rastair_env",
    "rasusa_env",
    "raven_env",
    "rclone_env",
    "regtools_env",
    "repeatmasker_env",
    "repeatmodeler_env",
    "resolvetaxonomy_env",
    "rgi_env",
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
    "samblaster_env",
    "samshee_env",
    "samtools_env",
    "sawfish_env",
    "scanpy_env",
    "scimap_env",
    "seacr_env",
    "sentieon_env",
    "seq2hla_env",
    "seqcluster_env",
    "seqfu_env",
    "seqkit_env",
    "seqtk_env",
    "sequali_env",
    "seqwish_env",
    "sgdemux_env",
    "shapeit5_env",
    "shinyngs_env",
    "sigprofiler_env",
    "simpleaf_env",
    "smoothxg_env",
    "sniffles_env",
    "snpeff_env",
    "sortmerna_env",
    "souporcell_env",
    "sourcepredict_env",
    "sourmash_env",
    "spaceranger_env",
    "spades_env",
    "sparsesignatures_env",
    "sratools_env",
    "star_env",
    "stardist_env",
    "starfusion_env",
    "stitch_env",
    "stranger_env",
    "strelka_env",
    "stringtie_env",
    "subread_env",
    "summarizedexperiment_env",
    "suppa_env",
    "survivor_env",
    "svanalyzer_env",
    "svim_env",
    "svtk_env",
    "svync_env",
    "sylph_env",
    "sylphtax_env",
    "taxonkit_env",
    "taxpasta_env",
    "tcoffee_env",
    "thermorawfileparser_env",
    "tiara_env",
    "tidk_env",
    "tinc_env",
    "toulligqc_env",
    "traitar_env",
    "transdecoder_env",
    "trgt_env",
    "trimgalore_env",
    "trimmomatic_env",
    "trinity_env",
    "trnascanse_env",
    "trust4_env",
    "truvari_env",
    "tx2gene_env",
    "tximport_env",
    "ucsc_env",
    "ultra_env",
    "umicollapse_env",
    "umitools_env",
    "unicycler_env",
    "untar_env",
    "untarfiles_env",
    "unzip_env",
    "unzipfiles_env",
    "upp_env",
    "variancepartition_env",
    "variantextractor_env",
    "vcf2maf_env",
    "vcflib_env",
    "vcftools_env",
    "viber_env",
    "vireo_env",
    "vrhyme_env",
    "vsearch_env",
    "vt_env",
    "wfmash_env",
    "wget_env",
    "wipertools_env",
    "wittyer_env",
    "xeniumranger_env",
    "xz_env",
    "yahs_env",
    "yara_env",
    "zip_env",
]
