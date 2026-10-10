# flyte-bio

Bioinformatics tools and pipelines for [Flyte](https://www.union.ai/docs/v2/flyte/).

> [!WARNING]
> This package is still experimental and is not ready for general use or
> publication yet. The APIs, module layout, and packaging details may still
> change.

flyte-bio has two layers:

- **Modules** (`flyte_bio.modules`) wrap individual tools as Flyte tasks:
  samtools, BWA, GATK, STAR, salmon, MultiQC and others. Each tool runs in its
  official biocontainer, so your own image doesn't need the tools installed.
- **Pipelines** (`flyte_bio.pipelines`) compose those tasks into complete
  workflows, ported from established Nextflow pipelines: RNA-seq
  quantification and germline variant calling.

## Install

flyte-bio isn't on PyPI yet. Install it from GitHub:

```bash
uv add "flyte-bio @ https://github.com/unionai-oss/flyte-bio/archive/refs/heads/main.zip"
```

It needs `flyte>=2.10.0`.

## Usage

Your code runs in your own Flyte task, which calls the flyte-bio tasks. That
needs two things:

- **flyte-bio installed in your task's image**, because your task imports it.
  The tool containers don't need it.
- **A `depends_on` on the flyte-bio environments**, so running or deploying
  your task also registers every tool task and its container.
  `flyte_bio.env` covers all of them.

```python
import flyte

from flyte_bio import env as bio_env

FLYTE_BIO = "flyte-bio @ https://github.com/unionai-oss/flyte-bio/archive/refs/heads/main.zip"

env = flyte.TaskEnvironment(
    name="my_analysis",
    image=flyte.Image.from_debian_base().with_pip_packages(FLYTE_BIO),
    resources=flyte.Resources(cpu=2, memory="4Gi"),
    depends_on=[bio_env],
)
```

### A module

Module functions are ordinary async calls inside your task:

```python
from flyte.io import File

from flyte_bio.modules.bedtools import bedtools_intersect


@env.task
async def overlaps(annotation: File, peaks: list[File]) -> File:
    return await bedtools_intersect(a=annotation, b=peaks, wa=True, f=0.5)
```

### A pipeline

A pipeline takes a samplesheet CSV, like the Nextflow pipeline's `--input`,
and can publish its results in that pipeline's `--outdir` layout:

```python
from flyte.io import Dir, File

from flyte_bio.pipelines.rnaseq import rnaseq


@env.task
async def quantify(samplesheet: File, fasta: File, gtf: File, outdir: str | None = None) -> Dir:
    result = await rnaseq(samplesheet, fasta=fasta, gtf=gtf, outdir=outdir, publish_results=True)
    return result.outdir
```

```bash
flyte run my_analysis.py quantify --samplesheet samples.csv --fasta genome.fa --gtf genes.gtf
```

[`examples/`](examples/) has a complete, runnable script for each pipeline.

## Documentation

- **[Modules](docs/modules.md):** how modules work, and every available tool.
- **[Pipelines](docs/pipelines/README.md):** inputs, outputs and the shared
  utilities, with a page per pipeline:
  - [rnaseq](docs/pipelines/rnaseq.md): STAR + salmon RNA-seq quantification.
  - [variant_calling](docs/pipelines/variant_calling.md): germline variant
    calling with BWA-MEM and GATK HaplotypeCaller.
- **[Contributing](docs/contributing.md):** porting a tool, and writing and
  running tests.
- **[Known issues](docs/known-issues.md):** temporary workarounds for
  upstream bugs.
