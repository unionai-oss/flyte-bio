# variant_calling

`flyte_bio.pipelines.variant_calling`: germline short-variant calling from
FASTQ, with BWA-MEM alignment, GATK preprocessing and GATK HaplotypeCaller.

It ports the germline HaplotypeCaller route of
[nf-core/sarek 3.10.0](https://nf-co.re/sarek/3.10.0/): upstream's default
preprocessing with `--tools haplotypecaller`.

## Example

```python
from flyte.io import Dir, File

from flyte_bio.pipelines.variant_calling import variant_calling


@env.task  # with flyte-bio in its image and depends_on=[flyte_bio.env]
async def call(
    samplesheet: File, fasta: File, dbsnp: File, known_indels: File, outdir: str | None = None
) -> Dir:
    result = await variant_calling(
        samplesheet, fasta=fasta, dbsnp=dbsnp, known_indels=[known_indels],
        outdir=outdir, publish_results=True,
    )
    return result.outdir
```

[`examples/variant_calling.py`](../../examples/variant_calling.py) is a
runnable version. To try it on the upstream test profile (one sample, two
lanes, a 40 kb slice of chr22):

```bash
DATA=https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome
curl -sO https://raw.githubusercontent.com/nf-core/sarek/3.10.0/tests/csv/3.0/fastq_single.csv
flyte run examples/variant_calling.py variant_calling_example \
    --samplesheet fastq_single.csv \
    --fasta $DATA/genome.fasta \
    --dbsnp $DATA/vcf/dbsnp_146.hg38.vcf.gz \
    --known_indels $DATA/vcf/mills_and_1000G.indels.vcf.gz \
    --intervals $DATA/genome.interval_list
```

The test reads come from the whole genome, so on that tiny reference they
produce no calls after base-quality recalibration (see [tests](#tests)). Real
data does.

## Samplesheet

The upstream columns, one row per lane:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
P1,XX,0,P1_N,L001,s3://my-bucket/P1_N_L001_R1.fastq.gz,s3://my-bucket/P1_N_L001_R2.fastq.gz
P1,XX,0,P1_N,L002,s3://my-bucket/P1_N_L002_R1.fastq.gz,s3://my-bucket/P1_N_L002_R2.fastq.gz
```

- **Grouping:** rows that share `patient` and `sample` are one sample.
  - Its rows must agree on `sex` and `status`.
  - Its lanes must be distinct.
- **`status`:** 0 is normal. Tumour samples (`status` 1) aren't supported yet.
- **`fastq_2`:** empty for single-end reads.
- **No whitespace** in patient, sample or lane IDs.

## What it runs

1. **Reference prep** (`prepare_reference`):
   - builds whatever isn't supplied: the FASTA index, sequence dictionary,
     BWA index, and the bgzip + tabix indexes of the known sites;
   - prepares the calling intervals: built from the FASTA index (one per
     sequence) when none are given, or an interval list converted to BED;
   - splits the intervals into scatter-gather chunks, sized by
     `nucleotides_per_second`.
2. **Mapping** (`map_reads`):
   - FastQC per lane;
   - BWA-MEM per lane with upstream's arguments (`-K 100000000 -Y`) and
     read group: `ID` is `<flowcell>.<sample>.<lane>`, with the flowcell taken
     from the FASTQ header when it has one. The output is coordinate-sorted.
3. **Duplicate marking** (`mark_duplicates`): GATK MarkDuplicates across the
   sample's lanes into `<sample>.md.cram`, then samtools stats and mosdepth
   (500 bp windows for WGS, the intervals with `wes`). With
   `skip_markduplicates`, the lanes are only merged, into
   `<sample>.sorted.cram`.
4. **Base-quality recalibration** (`recalibrate`):
   - BaseRecalibrator per interval chunk against dbSNP and the known indels,
     and the tables are gathered;
   - ApplyBQSR per chunk, merged into `<sample>.recal.cram`;
   - QC as for the duplicate-marked CRAM.
5. **Calling** (`call_variants`):
   - HaplotypeCaller per interval chunk (`--pcr-indel-model CONSERVATIVE`,
     annotating dbSNP IDs), and the chunks are merged into
     `<sample>.haplotypecaller.vcf.gz`;
   - then CNNScoreVariants and FilterVariantTranches against the known
     sites, giving `<sample>.haplotypecaller.filtered.vcf.gz`.
6. **VCF QC** (`vcf_qc`): bcftools stats, and vcftools Ts/Tv by count and by
   quality, and the FILTER summary.
7. **MultiQC** over every sample, with upstream's config (`multiqc_report`).

## Options

| Option | Default | Meaning |
|---|---|---|
| `fai`, `sequence_dict`, `bwa_index` | built | Prebuilt reference indexes |
| `dbsnp`, `dbsnp_tbi` | — | dbSNP (indexed if no `.tbi`). Recalibration and the filter need dbSNP and/or known indels |
| `known_indels`, `known_indels_tbi` | — | Known-indel VCFs (a list) |
| `intervals` | whole genome | BED or interval list of regions to call |
| `no_intervals` | `False` | Run every step over the whole genome in one go |
| `nucleotides_per_second` | 200000 | Scatter-chunk sizing: length / this ≈ runtime |
| `wes` | `False` | Exome/targeted data: mosdepth reports on the intervals |
| `seq_platform`, `seq_center` | `ILLUMINA`, — | Read-group `PL` and `CN` |
| `skip_markduplicates` | `False` | Merge lanes without marking duplicates |
| `skip_baserecalibrator` | `False` | Call on the duplicate-marked CRAM |
| `skip_haplotypecaller_filter` | `False` | Keep the unfiltered calls |
| `haplotypecaller_filter_args` | `--info-key CNN_1D` | FilterVariantTranches options |
| `gatk_pcr_indel_model` | `CONSERVATIVE` | HaplotypeCaller's PCR indel model |
| `skip_fastqc`, `skip_vcf_qc`, `skip_multiqc` | `False` | Turn QC steps off |
| `outdir`, `publish_results` | — | Publish the results ([outputs](README.md#outputs)) |

## Outputs

`variant_calling()` returns a `VariantCallingResult`:

- `reference`: the prepared reference and intervals.
- `samples`: per sample:
  - the mapped lanes and FastQC reports;
  - the duplicate-marked (or merged) CRAM with its QC;
  - the recalibrated CRAM, table and QC;
  - the calls: raw, CNN-scored and filtered, with `.final` being the one
    downstream steps use;
  - the VCF QC.
- `multiqc`.
- `outdir`: the published tree, when asked for.

Published results follow upstream's layout:

```
preprocessing/
  markduplicates/<sample>/<sample>.md.cram(.crai)   (mapped/ with skip_markduplicates)
  recal_table/<sample>/<sample>.recal.table
  recalibrated/<sample>/<sample>.recal.cram(.crai)
variant_calling/haplotypecaller/<sample>/
  <sample>.haplotypecaller.vcf.gz(.tbi)
  <sample>.haplotypecaller.filtered.vcf.gz(.tbi)
reports/
  fastqc/<sample>-<lane>/  markduplicates/<sample>/  samtools/<sample>/  mosdepth/<sample>/
  bcftools/haplotypecaller/<sample>/  vcftools/haplotypecaller/<sample>/
multiqc/
```

## Differences from upstream

- **Not ported:**
  - fastp trimming, and its default splitting of FASTQs into 50M-read
    chunks. Splitting changes how mapping is parallelised, not the
    alignments.
  - UMIs, and the bwa-mem2 and DRAGMAP aligners.
  - Joint germline calling, tumour/normal and tumour-only calling, and the
    other callers (DeepVariant, FreeBayes, Strelka, Manta, CNV callers and
    others).
  - Annotation (VEP, snpEff) and VCF post-processing.
  - Spark, Sentieon and Parabricks.
  - Starting from a later step (`--step`) with BAM or CRAM input.
- **A sample with no calls skips the CNN filter**, because GATK refuses to
  score or filter an empty VCF. Upstream fails the run.
- **The read group's `DS:` field** records the FASTA's storage path, where
  upstream records `--fasta` as given.
- **MultiQC:** the custom logo is dropped, and plot export is off for now
  ([known issues](../known-issues.md)).

## Tests

[`tests/pipelines/test_variant_calling.py`](../../tests/pipelines/test_variant_calling.py)
runs on the upstream test genome:
- **reference prep:** the built `.fai` matches upstream's byte for byte;
- **mapping:** upstream's read groups;
- **duplicate marking:** the test sample's two lanes are identical, so at
  least 50% of reads are marked duplicate;
- **BQSR,** scattered over two interval chunks;
- **calling:** upstream's HaplotypeCaller test from its mapped BAM, with the
  filter;
- **end to end:** the full `variant_calling()` route, published to Flyte
  storage and to an explicit `outdir`.

Two cases hit the tiny test data's limits, as upstream's own runs would:

- **No calls after recalibration:** the FASTQ test reads come from the whole
  genome, so on the 40 kb reference the reads that map carry about 8%
  mismatches. BQSR rightly lowers their qualities, and HaplotypeCaller calls
  nothing.
- **No filtering without recalibration:** the calls include indels that no
  known-indel site overlaps, which the tranche filter refuses.

The filter is exercised on upstream's mapped-BAM test instead.
