# rnaseq

`flyte_bio.pipelines.rnaseq`: RNA-seq read QC, STAR alignment, salmon
quantification, and alignment and expression QC.

It ports the default `star_salmon` route of
[nf-core/rnaseq 3.26.0](https://nf-co.re/rnaseq/3.26.0/), and its optional
UMI deduplication.

## Example

```python
from flyte.io import Dir, File

from flyte_bio.pipelines.rnaseq import rnaseq


@env.task  # with flyte-bio in its image and depends_on=[flyte_bio.env]
async def quantify(samplesheet: File, fasta: File, gtf: File, outdir: str | None = None) -> Dir:
    result = await rnaseq(samplesheet, fasta=fasta, gtf=gtf, outdir=outdir, publish_results=True)
    return result.outdir
```

[`examples/rnaseq.py`](../../examples/rnaseq.py) is a runnable version. To
try it on the upstream test profile (5 small yeast samples):

```bash
DATA=https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a
curl -sO $DATA/samplesheet/v3.10/samplesheet_test.csv
flyte run examples/rnaseq.py rnaseq_example \
    --samplesheet samplesheet_test.csv \
    --fasta $DATA/reference/genome.fasta \
    --gtf $DATA/reference/genes_with_empty_tid.gtf.gz \
    --transcript_fasta $DATA/reference/transcriptome.fasta
    # add --outdir s3://my-bucket/rnaseq to publish the results there
```

## Samplesheet

The upstream columns, one row per sequencing run:

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,s3://my-bucket/reads/ctrl1_R1.fastq.gz,s3://my-bucket/reads/ctrl1_R2.fastq.gz,auto
CONTROL_REP1,s3://my-bucket/reads/ctrl1b_R1.fastq.gz,s3://my-bucket/reads/ctrl1b_R2.fastq.gz,auto
TREATED_REP1,s3://my-bucket/reads/treat1_R1.fastq.gz,,reverse
```

- Rows that share a `sample` are that sample's runs. They're merged in file
  order, and they must agree on strandedness.
- `fastq_2` is empty for single-end reads. The column must be present anyway.
- `strandedness` is `forward`, `reverse`, `unstranded` or `auto`. With
  `auto`, salmon infers it from a 1M-read subsample, as upstream does.

## What it runs

1. **Reference prep** (`prepare_genome`):
   - decompress the references and filter the GTF to the genome's sequences;
   - append `additional_fasta` (e.g. spike-ins), if given;
   - derive the gene BED and the transcript FASTA, and build the STAR index.
   
   A GFF is converted to GTF first.
2. **Read preprocessing** (`preprocess_reads`):
   - merge a sample's runs, then `fq lint` the reads;
   - FastQC, then Trim Galore with FastQC on the trimmed reads;
   - drop samples left with fewer than `min_trimmed_reads`;
   - optionally BBSplit away reads from other genomes, then lint again.
3. **Strandedness inference** for `auto` samples (`infer_strandedness`).
4. **Alignment** with STAR, emitting a transcriptome BAM (`align_star`). The
   genome BAM is then sorted, indexed and summarised with samtools.
5. **Mapping gate:** samples below `min_mapped_reads` % uniquely mapped are
   still quantified, but skip every genome-BAM step.
6. **Duplicate marking** with Picard (`mark_duplicates`). With `with_umi`,
   UMI-tools deduplicates instead (`dedup_umi`).
7. **Alignment QC:**
   - dupRadar;
   - Qualimap, on a name-sorted copy;
   - RSeQC (7 modules);
   - featureCounts biotype QC;
   - StringTie;
   - bigWig coverage, per strand for stranded libraries.
8. **Quantification:** salmon in alignment mode on the transcriptome BAM
   (`quantify_salmon_bam`). Then tx2gene and tximport merge the samples into
   gene/transcript count, TPM and length matrices, and SummarizedExperiment
   RDS files (`merge_quantifications`).
9. **DESeq2 QC:** PCA and sample distances on the merged counts.
10. **MultiQC** with upstream's config and custom content
    (`multiqc_report`). The custom content covers failed samples, strandedness
    checks, merging paired-end samples, and FASTQ-to-sample name
    replacements.

## Options

`rnaseq()` takes the upstream parameters as keyword arguments:

| Option | Default | Meaning |
|---|---|---|
| `gtf` / `gff` | — | Gene annotation (one of them) |
| `transcript_fasta` | built with gffread | Transcript sequences |
| `additional_fasta` | — | Extra sequences appended to the genome and GTF |
| `salmon_index_dir` | built if needed | Salmon index, used for strandedness inference |
| `skip_bbsplit`, `bbsplit_fasta_list`, `bbsplit_index_dir` | skipped | BBSplit contamination removal: name → FASTA, or a prebuilt index |
| `skip_trimming`, `min_trimmed_reads` | trim, 10000 | Trim Galore, and the minimum reads after trimming |
| `skip_linting`, `extra_fqlint_args` | lint, `--disable-validator P001` | `fq lint` at each stage |
| `stranded_threshold`, `unstranded_threshold` | 0.8, 0.1 | Strandedness calls |
| `min_mapped_reads` | 5.0 | % uniquely mapped for the genome-BAM steps |
| `skip_markduplicates`, `skip_dupradar`, `skip_qualimap`, `skip_rseqc`, `rseqc_modules`, `skip_biotype_qc`, `skip_stringtie`, `skip_bigwig`, `skip_deseq2_qc`, `skip_fastqc`, `skip_multiqc` | all run | Turn individual steps off |
| `featurecounts_group_type`, `featurecounts_feature_type` | `gene_biotype`, `exon` | Biotype QC grouping |
| `with_umi` and the `umitools_*` / `umi_discard_read` options | off | UMI extraction and deduplication (below) |
| `seq_platform`, `seq_center` | — | Read-group fields in the STAR BAM |
| `outdir`, `publish_results` | — | Publish the results ([outputs](README.md#outputs)) |

### UMIs

`with_umi=True` enables upstream's UMI-tools route:
- **Extraction:** before trimming, UMIs move into read names, using
  `umitools_bc_pattern` and `umitools_extract_method`. Skip this with
  `skip_umi_extract` if your read names already carry UMIs.
- **Discarding a mate:** `umi_discard_read` (1 or 2) drops a mate that held
  only the UMI; the sample then continues as single-end.
- **Deduplication:** the genome BAM is deduplicated in place of
  MarkDuplicates. So is the transcriptome BAM, which is coordinate-sorted,
  deduplicated, name-sorted, and run through prepare-for-rsem for paired-end
  data, before salmon.

UMI patterns are passed as unquoted options, so they can't contain
whitespace.

## Outputs

`rnaseq()` returns an `RnaseqResult`:

- `genome`: the prepared references and STAR index.
- `samples`: per sample:
  - the preprocessing results and alignment;
  - salmon results, duplicate marking or UMI deduplication;
  - Qualimap, RSeQC, biotype QC, StringTie and bigWigs;
  - strandedness, and whether the mapping gate passed.
- `salmon`: the merged count/TPM/length matrices and SummarizedExperiment RDS
  files.
- `failed_trimming`: samples dropped by the trimming filter.
- `deseq2_qc`, `multiqc`.
- `outdir`: the published tree, when asked for.

Published results follow upstream's layout:

```
fastqc/{raw,trim}/  trimgalore/  fq_lint/  bbsplit/  umitools/
star_salmon/
  <sample>.markdup.sorted.bam(.bai)   (or <sample>.umi_dedup.sorted.bam)
  <sample>/                           salmon quantification
  salmon.merged.*                     gene/transcript matrices, tx2gene, .rds
  log/  samtools_stats/  picard_metrics/  qualimap/  rseqc/  dupradar/
  featurecounts/  stringtie/  bigwig/  deseq2_qc/  umitools/
multiqc/star_salmon/
```

## Differences from upstream

- **Not ported:**
  - the other aligner and quantifier routes (`star_rsem`, HISAT2, Bowtie2,
    pseudo-alignment with salmon or kallisto);
  - rRNA removal (SortMeRNA / RiboDetector), UMICollapse, Preseq, and
    Kraken2 / Bracken / Sylph;
  - fastp trimming, RustQC, Sentieon and Parabricks.
- **Reference prep:** without `transcript_fasta`, transcripts are extracted
  with gffread, where upstream uses RSEM's extractor.
- **MultiQC:** the Nextflow-specific report sections (run parameters,
  software versions, methods) are left out. The report includes samtools
  stats for both the sorted and the duplicate-marked BAM, where upstream
  includes one set. Plot export is off for now ([known issues](../known-issues.md)).
- **Sample IDs** can't contain whitespace.

## Tests

[`tests/pipelines/test_rnaseq.py`](../../tests/pipelines/test_rnaseq.py) runs
two cases:
- upstream's test profile, end to end;
- upstream's `--umi_dedup_tool umitools` test.

They check every output, and the paths in the published layout. Module
outputs are checked against upstream's md5s in the module tests.
