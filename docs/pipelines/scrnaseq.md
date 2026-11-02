# scrnaseq

`flyte_bio.pipelines.scrnaseq`: droplet single-cell RNA-seq from FASTQ to
count matrices, with simpleaf (alevin-fry) quantification, QCatch cell QC,
CellBender background removal, and h5ad and RDS outputs.

It ports the default `simpleaf` route of
[nf-core/scrnaseq 4.2.0](https://nf-co.re/scrnaseq/4.2.0/).

## Example

```python
from flyte.io import Dir, File

from flyte_bio.pipelines.scrnaseq import scrnaseq


@env.task  # with flyte-bio in its image and depends_on=[flyte_bio.env]
async def quantify(samplesheet: File, fasta: File, gtf: File, outdir: str | None = None) -> Dir:
    result = await scrnaseq(
        samplesheet, fasta=fasta, gtf=gtf, protocol="10XV3", outdir=outdir, publish_results=True
    )
    return result.outdir
```

[`examples/scrnaseq.py`](../../examples/scrnaseq.py) is a runnable version.
To try it on the upstream test profile (two 10x v2 samples against mouse
chr19):

```bash
DATA=https://github.com/nf-core/test-datasets/raw/scrnaseq
curl -sLO $DATA/samplesheet-2-0.csv
flyte run examples/scrnaseq.py scrnaseq_example \
    --samplesheet samplesheet-2-0.csv \
    --fasta $DATA/reference/GRCm38.p6.genome.chr19.fa \
    --gtf $DATA/reference/gencode.vM19.annotation.chr19.gtf \
    --protocol 10XV2 \
    --skip_qcatch --skip_cellbender
```

QCatch and CellBender need a realistic number of cells, so upstream's test
profile skips them. Leave them on for real data.

## Samplesheet

The upstream columns, one row per lane:

```csv
sample,fastq_1,fastq_2,expected_cells,seq_center
S1,s3://my-bucket/S1_L001_R1.fastq.gz,s3://my-bucket/S1_L001_R2.fastq.gz,5000,Broad Institute
S1,s3://my-bucket/S1_L002_R1.fastq.gz,s3://my-bucket/S1_L002_R2.fastq.gz,5000,Broad Institute
```

- **`fastq_1`** carries the cell barcode and UMI, **`fastq_2`** the cDNA.
  Both are required.
- **`expected_cells`** and **`seq_center`** are optional. Rows of one sample
  must agree on them. Every column is copied onto the cells of the combined
  matrix.
- **Grouping:** rows that share `sample` are one sample's lanes.
- **No whitespace** in sample IDs.

## What it runs

1. **Reference prep** (`prepare_reference`): decompresses the FASTA and GTF
   if gzipped, and keeps the GTF lines on the genome's sequences
   (`<genome>_genes.gtf`).
2. **Index** (`prepare_index`): the simpleaf (piscem) index of the spliced +
   intronic transcriptome, with its transcript-to-gene table. A prebuilt
   `simpleaf_index` (with `txp2gene`) is used as given, and
   `transcript_fasta` + `txp2gene` build it from transcripts instead.
3. **Read QC** (`read_qc`): FastQC over all of a sample's FASTQs, as
   `<sample>_<n>_fastqc.*` in lane order (R1, R2, ...).
4. **Quantification** (`quantify`): simpleaf quant per sample with the
   protocol's chemistry, against its barcode whitelist (fetched from
   upstream's assets, or `barcode_whitelist`), unfiltered (`unfiltered-pl`),
   with `cr-like` UMI resolution.
5. **Cell QC** (`cell_qc`): QCatch cell calling and its HTML report,
   filtered h5ad and metrics summary.
6. **Matrices** (`convert_matrices`): each sample's raw (not cell-filtered)
   counts as `<sample>_raw_matrix.h5ad`, all of them combined into
   `combined_raw_matrix.h5ad` with the samplesheet columns, and each as
   Seurat and SingleCellExperiment RDS.
7. **Background removal** (`remove_background`): CellBender calls cells on
   each raw matrix, on the CPU or, with `cellbender_gpu`, a GPU. Each raw
   matrix restricted to its cells becomes
   `<sample>_cellbender_filter_matrix.h5ad`, combined and converted like the
   raw ones. Those keep the raw counts; CellBender's denoised counts are in
   its own `.h5` outputs.
8. **MultiQC** over the FastQC reports and simpleaf quantifications, with
   upstream's config (`multiqc_report`).

## Options

| Option | Default | Meaning |
|---|---|---|
| `protocol` | required | `10XV1`, `10XV2`, `10XV3`, `10XV4` or `dropseq`. Any other value goes to simpleaf as its chemistry, e.g. a custom read geometry |
| `barcode_whitelist` | the protocol's | Cell barcode whitelist |
| `simpleaf_index`, `txp2gene` | built | A prebuilt index and its transcript-to-gene table |
| `transcript_fasta` | — | Build the index from transcripts (needs `txp2gene`) |
| `simpleaf_umi_resolution` | `cr-like` | alevin-fry's UMI resolution |
| `skip_qcatch` | `False` | Skip QCatch |
| `remove_doublets`, `qcatch_n_partitions` | `False`, — | QCatch doublet removal; its partition count, when it can't infer the chemistry |
| `skip_cellbender` | `False` | Skip CellBender |
| `cellbender_gpu` | `False` | Run CellBender on a GPU (Broad's CUDA image, one GPU) |
| `cellbender_args` | — | Extra CellBender options, e.g. `--epochs 150` |
| `skip_fastqc`, `skip_multiqc` | `False` | Turn QC steps off |
| `save_reference` | `False` | Also publish the filtered GTF and the built index |
| `save_align_intermeds` | `True` | Publish simpleaf's mapping (`af_map`) |
| `outdir`, `publish_results` | — | Publish the results ([outputs](README.md#outputs)) |

## Outputs

`scrnaseq()` returns a `ScrnaseqResult`:

- `reference` and `index`: the prepared reference and simpleaf index.
- `samples`: per sample, the FastQC reports, the simpleaf mapping and
  quantification, and the QCatch results.
- `raw`: the raw matrices, per sample and combined.
- `cellbender`: per sample, CellBender's outputs, and its cells' matrices.
- `multiqc`.
- `outdir`: the published tree, when asked for.

Published results follow upstream's layout:

```
fastqc/<sample>_<n>_fastqc.{html,zip}
simpleaf/
  <sample>/
    simpleaf_quant/af_map/     (unless save_align_intermeds=False)
    simpleaf_quant/af_quant/   (the counts in alevin/)
    qcatch/
    cellbender_removebackground/
  mtx_conversions/
    <sample>/<sample>_{raw,cellbender_filter}_matrix.{h5ad,sce.rds,seurat.rds}
    combined_{raw,cellbender_filter}_matrix.{h5ad,sce.rds,seurat.rds}
  simpleaf_index/              (save_reference)
gtf_filter/                    (save_reference)
multiqc/
```

## Differences from upstream

- **Not ported:** the other aligners (STARsolo, kallisto | bustools, Cell
  Ranger, Cell Ranger ARC and multi), iGenomes references (`--genome`), and
  MultiQC's custom config, title and logo options.
- **`cellbender_args`** passes options to CellBender, which upstream does
  through `ext.args`.
- **MultiQC:**
  - it's 1.33, upstream's 1.34;
  - the report has no software-versions, workflow-summary or methods
    sections, which come from Nextflow;
  - of each `af_quant`, only its JSON summaries go to MultiQC. No MultiQC
    module reads alevin-fry's output, so the report is the same.
- **`pipeline_info/`** isn't published.
- **The h5ad and RDS files' bytes** differ from upstream's snapshot, because
  they record the tools' command lines and file paths. The tests check their
  content instead.

## Tests

[`tests/pipelines/test_scrnaseq.py`](../../tests/pipelines/test_scrnaseq.py)
runs on the upstream test profile:

- **reference prep, read QC and quantification**, checking the GTF filter,
  the FastQC report names and that both of a sample's lanes are mapped;
- **matrices:** each h5ad's shape and total count match simpleaf's own, with
  sorted barcodes, unversioned gene IDs and sample labels, and the combined
  matrix holds every sample's cells;
- **MultiQC:** every FastQC report is in the report, and the plots are
  exported;
- **CellBender,** on upstream's module test matrix (a raw PBMC sample) for 5
  epochs, since it can't train on the profile's few cells: the filtered h5ad
  holds exactly CellBender's cells, with their raw counts;
- **end to end:** the whole `scrnaseq()` route, published to Flyte storage,
  gives exactly the files of upstream's simpleaf test snapshot (outside
  MultiQC's internals); with FastQC and MultiQC skipped, the reference saved
  and the mapping left out, published to an explicit `outdir`.

QCatch fails on the profile's few cells, as upstream's test profile expects,
so its module test runs on QCatch's own 1k-PBMC data.
