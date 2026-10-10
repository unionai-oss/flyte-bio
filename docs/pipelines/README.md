# Pipelines

`flyte_bio.pipelines` composes [modules](../modules.md) into complete
workflows. Each one ports a published Nextflow pipeline's default route and
tests against that pipeline's own test profile.

| Pipeline | What it does | Ported from |
|---|---|---|
| [rnaseq](rnaseq.md) | RNA-seq: STAR alignment, salmon quantification, and alignment and expression QC | rnaseq 3.26.0, `star_salmon` route |
| [variant_calling](variant_calling.md) | Germline short variants: BWA-MEM, GATK preprocessing and HaplotypeCaller | sarek 3.10.0, germline HaplotypeCaller route |

## How pipelines run

A pipeline is a set of plain async Python functions, not tasks. You call them
from your own task, which then fans out to the module tasks:

```python
from flyte.io import Dir, File

from flyte_bio.pipelines.variant_calling import variant_calling


@env.task  # with flyte-bio in its image and depends_on=[flyte_bio.env]
async def call(samplesheet: File, fasta: File, dbsnp: File, known_indels: File) -> Dir:
    result = await variant_calling(
        samplesheet, fasta=fasta, dbsnp=dbsnp, known_indels=[known_indels], publish_results=True
    )
    return result.outdir
```

Your task does little work itself: it reads the samplesheet, decides what to
run and collects results. The default `cpu=2, memory="4Gi"` is plenty. Each
tool runs as its own task, in parallel across samples, with Flyte's caching:
re-running with the same inputs reuses every step that hasn't changed.

Each pipeline's main function (`rnaseq()`, `variant_calling()`) runs the whole
route. Its steps are public too, e.g. `prepare_reference`, `map_reads` and
`call_variants`, for when you need only part of it.

## Inputs: the samplesheet

Like the Nextflow pipelines' `--input`, a pipeline takes a CSV samplesheet
with one row per sequencing run or lane. The columns are the upstream
pipeline's, documented on each pipeline's page. The samplesheet is checked
before anything runs, so a missing column fails the run with an error naming
the column.

- **FASTQ paths** must be URIs the cluster can read: `s3://`, `gs://`,
  `https://`. They're used where they are, not copied.
- **The samplesheet itself** is read by your task, so give it as an
  `s3://`/`gs://` URI or a local file (`flyte run` uploads it). Reading
  `https://` files inside a task isn't reliable.

Reference inputs (FASTA, GTF, known sites) are task inputs like any other
`File`. Indexes that aren't supplied are built, and cached for the next run.

## Outputs

A pipeline returns a result object holding every output as Flyte `File`s and
`Dir`s, e.g. `RnaseqResult` or `VariantCallingResult`. Everything stays in
Flyte's storage.

Like the Nextflow pipelines' `--outdir`, a pipeline can also lay its results
out in the upstream folder structure:

- `outdir="s3://my-bucket/run-1"` publishes to your bucket. The cluster's
  role needs write access.
- `publish_results=True` publishes to Flyte's own storage.

Either way, `result.outdir` is that tree as a `Dir`, which a downstream task
can take as an input. Copies within one object store are server-side. As
upstream does by default, intermediates aren't published.

## Shared utilities

Two generic helpers are available for writing new pipelines:

- **`flyte_bio.samplesheet.read_samplesheet(file)`** returns a CSV's rows as
  dicts, in file order, with values stripped. A pipeline indexes the columns
  it defines; a missing column raises an error naming the column and line.
- **`flyte_bio.publish.publish(layout, outdir=None)`** copies a mapping of
  relative path → `File`/`Dir` into one tree and returns it as a `Dir`. It
  rejects two outputs mapped to the same path. Each pipeline's
  `*_layout(result)` function builds its mapping.

## Running the examples

[`examples/`](../../examples/) has a script per pipeline. Its docstring has a
command to run it on the upstream test data:

```bash
flyte run examples/rnaseq.py rnaseq_example --samplesheet ... --fasta ... --gtf ...
```

Run from the repository root, against whichever cluster your Flyte config
points at.
