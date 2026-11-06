# Module porting plan

The order in which we port the remaining upstream modules, most used first.

## How modules are ranked

A module's popularity is the number of released pipelines that use it. The
counts come from the module list of each pipeline's latest release in the
nf-core pipeline index (`https://nf-co.re/pipelines.json`), taken on
2026-11-05. That index covers 94 released pipelines that use shared modules.
Pipelines written before shared modules existed aren't counted.

[`module-popularity.csv`](module-popularity.csv) has every used module with its
count, tranche, whether it's ported, and the pipelines that use it. Its
`ported` column is kept current as modules land; the counts below are from
the ranking.

Not ranked:

- **Software-version reporting** (`custom/dumpsoftwareversions`): it's Nextflow
  bookkeeping, and Flyte records task images itself.
- **Pipelines' own local modules**, which aren't shared. Porting a whole
  pipeline also means porting its local modules, as the rnaseq, scrnaseq and
  variant_calling ports did.

## Where we are

Of 2,128 upstream modules, 1,043 are used by at least one released pipeline
(not counting version reporting). We've ported 71 of those, including 12 of
the 13 used by 20 or more pipelines.

| Tranche | Pipelines using each module | Modules | Already ported | To port | Families to port | Mean pipeline coverage after it |
|---|---|---|---|---|---|---|
| 1 | 10 or more | 29 | 21 | 8 | 7 | 42% |
| 2 | 5–9 | 49 | 13 | 36 | 24 | 51% |
| 3 | 3–4 | 86 | 18 | 68 | 42 | 60% |
| 4 | 2 | 152 | 9 | 143 | 95 | 70% |
| 5 | 1 | 727 | 10 | 717 | 371 | 100% |

A pipeline's coverage is the share of its shared modules that are ported: 37%
on average today. Today 4 pipelines have every shared module ported. That's 5 after tranche 1, 12 after
tranche 3 and 16 after tranche 4. The long tail is what completes a pipeline.

## Recommendation

1. **Tranches 1–3 by popularity:** 112 modules in 64 families. Each
   module here serves at least 3 pipelines, and together they bring the
   average pipeline from 37% covered today to 60%.
2. **After tranche 3, port whole pipelines.** Tranches 4 and 5 are 860 modules,
   most used by a single pipeline. Ported in popularity order, they would
   complete few pipelines: tranche 4's 143 modules complete only 4 more. Ported per target pipeline, each
   batch finishes one, local modules included. Pick target pipelines by demand.
3. **Port the 1,084 modules that no released pipeline uses only when
   someone needs one.**

Port a whole family at once, even across tranches. One PR per family, e.g.
bwamem2's `index` (tranche 1) with its `mem` (tranche 2), since they share an
image, an environment and test data.

## Tranche 1: 10 or more pipelines

| Family | Modules (pipelines) | Notes |
|---|---|---|
| fastp | fastp (21) | Trimming and QC, often instead of Trim Galore. JSON and HTML reports feed MultiQC |
| bowtie2 | build (18), align (17) | Short-read aligner. align pipes through samtools, like bwa mem |
| gawk | gawk (13) | A generic awk step whose program each pipeline supplies. Port as a shell task taking the program |
| kraken2 | kraken2 (11) | Taxonomic classification. Needs a database: test with upstream's small test database |
| minimap2 | align (11) | Long-read and assembly aligner. Take `index` (tranche 2) too |
| bwamem2 | index (11) | Take `mem` (tranche 2) too. Indexing is memory-heavy |
| cat | cat (10) | Concatenate files, decompressing as needed |

## Tranche 2: 5–9 pipelines

| Family | Modules (pipelines) |
|---|---|
| bcftools | sort (9), norm (9), view (8), filter (6), annotate (6), mpileup (5), merge (5) |
| tabix | tabix (9), bgzip (9) |
| bwamem2 | mem (9) |
| find | concatenate (9) |
| unzip | unzip (9) |
| cutadapt | cutadapt (7) |
| custom | getchromsizes (7) |
| picard | collectmultiplemetrics (7) |
| nanoplot | nanoplot (7) |
| prokka | prokka (7) |
| quast | quast (7) |
| samtools | fastq (7), collatefastq (5) |
| ensemblvep | vep (7), download (6) |
| hmmer | hmmsearch (6) |
| mafft | align (6) |
| pigz | compress (6), uncompress (5) |
| preseq | lcextrap (6) |
| busco | busco (6) |
| minimap2 | index (6) |
| spades | spades (6) |
| bakta | baktadbdownload (5), bakta (5) |
| gatk4 | bedtointervallist (5) |
| bbmap | bbduk (5) |
| seqkit | stats (5), replace (5) |

## Tranche 3: 3–4 pipelines

| Family | Modules (pipelines) |
|---|---|
| gatk4 | variantfiltration (4), mutect2 (4), mergemutectstats (3), intervallisttools (3), filtermutectcalls (3) |
| bcftools | index (4), concat (4), reheader (3), consensus (3) |
| picard | createsequencedictionary (4), addorreplacereadgroups (4), collecthsmetrics (3) |
| mmseqs | createtsv (4), createdb (4), linclust (3), cluster (3) |
| hisat2 | build (4), extractsplicesites (3), align (3) |
| hmmer | hmmrank (3), hmmbuild (3), hmmalign (3) |
| deeptools | plotheatmap (3), plotfingerprint (3), computematrix (3) |
| samtools | convert (4), dict (3) |
| bowtie | build (4), align (3) |
| porechop | porechop (4), abi (3) |
| csvtk | join (4), concat (3) |
| snpeff | snpeff (4), download (3) |
| dragmap | hashtable (4), align (3) |
| bedtools | maskfasta (3), getfasta (3) |
| blast | makeblastdb (3), blastn (3) |
| seqkit | seq (4) |
| others, one module each | duckdb/table2parquet (4), flye (4), megahit (4), trimmomatic (4), seqtk/sample (4), sortmerna (4), prodigal (4), deepvariant/rundeepvariant (4), stringtie/merge (4), rsem/preparereference (4), clustalo/align (3), vsearch/cluster (3), macs2/callpeak (3), filtlong (3), toulligqc (3), rseqc/tin (3), sourmash/sketch (3), krona/ktimporttaxonomy (3), ucsc/wigtobigwig (3), rtgtools/format (3), bbmap/bbnorm (3), cellpose (3), seqfu/stats (3), kallisto/index (3), umicollapse (3), bracken (3) |

## What makes a module harder

A rough guide for sizing and ordering work within a tranche:

- **Small utilities** (gawk, cat/cat, find/concatenate, unzip, pigz, csvtk,
  duckdb). Quick to port. A few might be simpler as one generic helper than as
  separate tools.
- **Reference databases** (kraken2, bracken, ensemblvep, snpeff, busco, bakta,
  blast, sortmerna). The tool needs a database as an input. Upstream's tests use
  small test databases, and we should too. Downloading real databases
  (`ensemblvep/download`, `snpeff/download`, `bakta/baktadbdownload`) should be
  tasks of their own, so the result is cached.
- **Heavy compute** (spades, megahit, flye, deepvariant, bwamem2 index,
  dragmap hashtable). Size the default resources from upstream's labels. The
  tool's threads and memory come from `flyte_bio.limits`, so an override
  reaches the tool. DeepVariant and Cellpose can use GPUs.
- **Long reads** (minimap2, nanoplot, porechop, flye, filtlong, toulligqc).
  New test data, from upstream's own test datasets.

## Done means

The same as for every module so far ([contributing](contributing.md#porting-a-module)):

1. a shell task in the container upstream pins;
2. tests against upstream's own test cases, with md5s wherever the output is
   reproducible;
3. a row in the [modules table](modules.md#available-modules).

## Refreshing the ranking

Re-rank before each tranche, since pipelines change their modules between
releases. Count, for each module, the released pipelines whose latest release
lists it in `components.modules` in `pipelines.json`. A module listed as
`samtools_sort` is the module at `modules/nf-core/samtools/sort` in the modules
repository.
