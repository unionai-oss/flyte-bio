# Pipelines

`flyte_bio.pipelines` composes [modules](../modules.md) into complete
workflows. Each one ports a published Nextflow pipeline's default route and
tests against that pipeline's own test profile.

| Pipeline | What it does | Ported from |
|---|---|---|
| [rnaseq](rnaseq.md) | RNA-seq: STAR alignment, salmon quantification, and alignment and expression QC | rnaseq 3.26.0, `star_salmon` route |
| [scrnaseq](scrnaseq.md) | Droplet single-cell RNA-seq: simpleaf (alevin-fry) quantification, cell QC, CellBender, h5ad/RDS matrices | scrnaseq 4.2.0, `simpleaf` route |
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

Each pipeline's main function (`rnaseq()`, `scrnaseq()`, `variant_calling()`) runs the whole
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
`Dir`s, e.g. `RnaseqResult`, `ScrnaseqResult` or `VariantCallingResult`. Everything stays in
Flyte's storage.

Like the Nextflow pipelines' `--outdir`, a pipeline can also lay its results
out in the upstream folder structure:

- `outdir="s3://my-bucket/run-1"` publishes to your bucket. The cluster's
  role needs write access.
- `publish_results=True` publishes to Flyte's own storage.

Either way, `result.outdir` is that tree as a `Dir`, which a downstream task
can take as an input. Copies within one object store are server-side. As
upstream does by default, intermediates aren't published.

## Configuring a run with a file

Instead of passing options as arguments, you can keep them in a YAML file,
like a Nextflow `-params-file`. `flyte_bio.config` turns a config into a
pipeline's keyword arguments. It needs `flyte-bio[config]` in your task's
image.

A config is flat: one key per option of the pipeline function, with the same
name and default. Files and directories are URIs:

```yaml
samplesheet: s3://my-bucket/samples.csv
protocol: 10XV3
fasta: s3://my-bucket/ref/genome.fa
gtf: s3://my-bucket/ref/genes.gtf
skip_cellbender: true
```

Your task takes the config as an OmegaConf `DictConfig` and calls
`run_pipeline`:

```python
from flyte.io import Dir, File
from omegaconf import DictConfig

from flyte_bio.config import run_pipeline
from flyte_bio.pipelines.scrnaseq import scrnaseq


@env.task
async def run(cfg: DictConfig, samplesheet: File | None = None) -> Dir:
    result = await run_pipeline(scrnaseq, cfg, samplesheet=samplesheet, publish_results=True)
    return result.outdir
```

The config is checked against the function's signature before anything runs.
A misspelled option, a value of the wrong type or a missing required option
fails with an error naming it. Keyword arguments to `run_pipeline` take
precedence over the config; `None` ones are ignored. That's how a local file
reaches the pipeline: `flyte run` uploads a task's `File` arguments, but a
local path in a config is just a string, so the cluster can't read it.

### With Hydra

[`src/flyte_bio/conf/`](../../src/flyte_bio/conf/) has a
[Hydra](https://hydra.cc) config per pipeline. `<pipeline>/<pipeline>.yaml`
lists every option with its default, with `???` for the required ones, and
`profile=test` adds upstream's test data. `flyte hydra run` composes a config
and passes it to your task. It needs `flyte-bio[hydra]` where you launch the
run:

```bash
flyte hydra run --config-path src/flyte_bio/conf/scrnaseq --config-name scrnaseq \
    examples/scrnaseq.py scrnaseq_from_config \
    --samplesheet samplesheet.csv --cfg profile=test --cfg skip_fastqc=true
```

`--cfg key=value` sets any option. With `--multirun`, comma-separated values
launch one run per value (`--cfg simpleaf_umi_resolution=cr-like,parsimony`).

Your own config can build on flyte-bio's, so it only lists what differs. Put
this in `conf/my_scrnaseq.yaml` and pass `--config-path conf --config-name
my_scrnaseq`:

```yaml
hydra:
  searchpath:
    - pkg://flyte_bio.conf
defaults:
  - scrnaseq/scrnaseq@_here_
  - _self_
protocol: 10XV3
fasta: s3://my-bucket/ref/genome.fa
gtf: s3://my-bucket/ref/genes.gtf
```

Options that aren't arguments of the pipeline function aren't configurable
yet. In particular, per-tool CPU and memory aren't, because the module
functions don't take resources.

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
