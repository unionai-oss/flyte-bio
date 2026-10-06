# flyte-bio

Curated bioinformatics building blocks for Flyte.

> [!WARNING]
> This package is still experimental and is not ready for general use or
> publication yet. The APIs, module layout, and packaging details may still
> change.

A collection of typed Flyte tasks for popular bioinformatics tools — bedtools,
samtools, bcftools, GATK, Picard, BWA, STAR, salmon, kallisto, and friends —
each shipped via its official biocontainer image, plus higher-level pipelines
composed from those tasks. Built on `flyte.extras.shell`.

The package is split into two sub-packages:

- `flyte_bio.modules` — one file per tool family. Each exposes callable shell
  tasks (`bedtools_intersect`, `samtools_view`, ...) and a module-level `env`
  (a `flyte.TaskEnvironment`) that pipelines add to their `depends_on`. Adding
  the env causes the deploy pipeline to register the task and pull/build the
  underlying biocontainer image.
- `flyte_bio.pipelines` — higher-level workflows (`rnaseq`, variant calling,
  ...) that compose module tasks. Each pipeline exposes its own
  `TaskEnvironment`.

The package root re-exports `flyte_bio.env`, an aggregate environment that
depends on every module env. Pipelines that want broad access to the bio
module suite can depend on this one env.

## Install

This package is not yet published. For development, clone and install with
[uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/<org>/flyte-bio
cd flyte-bio
uv sync
```

`flyte` is currently pinned to a dev branch of
[flyteorg/flyte-sdk](https://github.com/flyteorg/flyte-sdk) (PR
[#1055](https://github.com/flyteorg/flyte-sdk/pull/1055)), which adds the
`flyte.extras.shell` task wrapper this package builds on.

## Usage

```python
import flyte
from flyte.io import File
from flyte_bio import env as bio_env
from flyte_bio.modules.bedtools import bedtools_intersect

env = flyte.TaskEnvironment(
    name="genomics_pipeline",
    depends_on=[bio_env],
)

@env.task
async def pipeline(annotation: File, peaks: list[File]) -> list[File]:
    return await bedtools_intersect(a=annotation, b=peaks, wa=True, f=0.5)
```

## Currently available

### Modules (`flyte_bio.modules`)

- `bbmap` — `bbsplit_index`, `bbsplit`
- `bedtools` — `bedtools_intersect`, `bedtools_sort`, `bedtools_merge`
- `cat` — `cat_fastq`
- `catadditionalfasta` — `cat_additional_fasta`
- `dupradar` — `dupradar`
- `fastqc` — `fastqc`
- `fq` — `fq_subsample`, `fq_lint`
- `gffread` — `gffread_gff_to_gtf`, `gffread_transcripts_fasta`
- `gtf2bed` — `gtf2bed`
- `gtffilter` — `gtf_filter`
- `gunzip` — `gunzip`
- `picard` — `picard_markduplicates`
- `qualimap` — `qualimap_rnaseq`
- `salmon` — `salmon_index`, `salmon_quant_reads`, `salmon_quant_bam`
- `samtools` — `samtools_faidx`, `samtools_sort`, `samtools_index`, `samtools_stats`, `samtools_flagstat`, `samtools_idxstats`
- `star` — `star_genome_generate`, `star_align`
- `summarizedexperiment` — `summarized_experiment`
- `trimgalore` — `trimgalore`
- `tx2gene` — `tx2gene`
- `tximport` — `tximport`, `collect_quants`
- `untar` — `untar`

### Pipelines (`flyte_bio.pipelines`)

- `rnaseq` — STAR alignment + salmon quantification (the `star_salmon` path
  of rnaseq 3.26.0): `prepare_genome`, `align_star`, `quantify_salmon_bam`,
  `rnaseq`, with read linting, QC/trimming, BBSplit, strandedness inference,
  duplicate marking, dupRadar, Qualimap and the tximport/SummarizedExperiment
  merge. StringTie, bigWig coverage, RSeQC, biotype QC, deseq2_qc and MultiQC
  are not ported yet.

More tools and pipelines are added as needed. Contributions following the same
pattern (one file per tool family, sharing one biocontainer image, exposing a
module-level `env`) are welcome once the plugin is ready to stabilize. The
top-level `flyte_bio.env` is intended to grow alongside the modules.

## Workarounds pending upstream fixes

### Optional File inputs are `list[File]` (flyteorg/flyte#8118)

Optional file inputs that should be `File | None` are declared as `list[File]`
(0 or 1 item, default `[]`). On the cluster, copilot stages a *set*
`Optional[File]` as a bare path instead of the per-input directory that shell
scripts glob, so the input was silently ignored. Revert once
[flyteorg/flyte#8118](https://github.com/flyteorg/flyte/pull/8118) is
**merged and deployed** (for Union clusters: the `flyte2` submodule pin in the
`cloud` repo includes it and the cluster has been redeployed). Merging alone
doesn't change what the cluster runs.

To revert:

1. Inputs back to `File | None` and drop their `[]` defaults:
   `salmon_index.genome_fasta`, `salmon_quant_reads.reads_2`,
   `star_align.reads_2`, `bedtools_intersect.g`, `bedtools_sort.g`,
   `summarized_experiment_cmd.rowdata` / `.coldata`, `fastqc.reads_2`,
   `trimgalore_cmd.reads_2`, `fq_subsample_cmd.reads_2`, `fq_lint.reads_2`,
   `bbsplit_cmd.reads_2`, `picard_markduplicates_cmd.fasta` / `.fai`
   (each is marked with a `flyteorg/flyte#8118` comment).
2. Update the scripts that read them through a `nullglob` array
   (salmon index/quant, STAR, summarizedexperiment, fastqc, trimgalore, fq,
   bbsplit, picard); the
   bedtools `-g` flag
   needs no script change.
3. Update callers: `align_star` in `pipelines/rnaseq.py` (`reads_2=[...]`),
   the `summarized_experiment`, `trimgalore`, `fq_subsample`, `bbsplit` and
   `picard_markduplicates` wrappers,
   the `lint` helper and the `fastqc` call in
   `pipelines/rnaseq.py`, `tests/modules/test_fastqc.py`,
   `tests/modules/test_star.py` and `tests/modules/test_salmon.py`.
4. Rerun the suite. `test_index` asserts the salmon index really has decoys
   and the rnaseq pipeline test exercises paired-end STAR, so a staging
   regression fails loudly.

## Porting an nf-core module

Most wrappers in this package are direct ports of [nf-core/modules](https://github.com/nf-core/modules):
same biocontainer, same flag set, same test fixtures. Following the same source
keeps semantics consistent and gives us a free regression oracle — every
nf-core module ships a snapshot file with expected output MD5s that we can
assert against.

The pattern, end to end:

### 1. Pull the upstream module

For a tool like `samtools sort` (`modules/nf-core/samtools/sort/`):

```bash
gh api 'repos/nf-core/modules/contents/modules/nf-core/samtools/sort/main.nf' \
  --jq '.content' | base64 -d
```

The three files worth reading:

- `main.nf` — gives the **pinned biocontainer URI** (the `quay.io/biocontainers/...`
  line) and the **invocation template** showing which flags are exposed.
- `tests/main.nf.test` — names the upstream test cases and which fixtures
  (from [nf-core/test-datasets@modules](https://github.com/nf-core/test-datasets/tree/modules))
  they consume.
- `tests/main.nf.test.snap` — holds the expected output MD5 per test case
  (look for `"<filename>:md5,<hex>"`). These are our assertion targets.

### 2. Write the shell wrapper

Add a new module under `src/flyte_bio/modules/`, following `bedtools.py` as
the template. Each tool family lives in one file and shares one biocontainer
image:

```python
from flyte.extras import shell
from flyte.io import File
import flyte

SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"  # from main.nf
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")  # nf-core process_single

samtools_sort = shell.create(
    name="samtools_sort",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "input": File,
        "n": bool,         # sort by read name
        "threads": int | None,
        ...
    },
    defaults={"n": False, ...},  # every non-optional bool needs a default
    outputs={"bam": File},
    script=r"""
        samtools sort {flags.n} {flags.threads} \
            -o {outputs.bam} {inputs.input}
    """,
)

env = flyte.TaskEnvironment.from_task("samtools", samtools_sort.as_task(), ...)
```

A few rules worth internalizing:

- **Image URI** comes verbatim from the upstream `main.nf` container line
  (the `quay.io/biocontainers/...` branch, not the singularity URL).
- **Default resources** mirror nf-core's `process_single` label (1 cpu, 6 GiB)
  unless the upstream module uses a different one (`process_medium`,
  `process_high`, ...). Override per-call when a workload needs more.
- **Defaults** must cover every bool input. Optional scalars (`T | None`)
  are implicitly `None` and don't need an entry. Without defaults, callers
  have to pass every flag explicitly.
- **Flag names** default to `-{input_name}`. Use `flag_aliases={"py_name": "-cli-name"}`
  when the CLI flag would collide with a Python keyword or an output name.
- **Bundle** every task in a module-level `env` via `flyte.TaskEnvironment.from_task`,
  and add that env to the aggregate `flyte_bio.modules.env` in
  `src/flyte_bio/modules/__init__.py`.

### 3. Write an end-to-end test

Tests in this repo are themselves Flyte tasks. Locally they run via the
local executor; remotely (`flyte run`) each test lands in its own pod, so
hundreds of tests scale horizontally across the cluster instead of
contending for one machine's docker daemon. There is no pytest — the
"runner" is just an `asyncio.gather` inside a top-level task.

The helpers in [`tests/framework.py`](tests/framework.py) cover everything
you'll need:

- `env` — the `TaskEnvironment` to decorate tests with. Already depends on
  `flyte_bio.modules.env`, so every wrapped tool is in scope.
- `fixture(rel_path)` — cached Flyte task that fetches a file from
  `nf-core/test-datasets@modules`. The path matches what an upstream
  `main.nf.test` references via `params.modules_testdata_base_path` —
  paste it straight in.
- `assert_md5(file, expected, label=...)` — streams the file and raises
  `AssertionError` on md5 mismatch.
- `gather_tests([...])` — runs tests in parallel, prints a summary,
  raises if any failed. Use it inside a driver task.

A new test is one async function:

```python
from flyte_bio.modules.samtools import samtools_sort
from tests.framework import env, fixture, assert_md5

@env.task
async def test_samtools_sort() -> None:
    # nf-core: samtools/sort "test_samtools_sort"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.bam")
    out = await samtools_sort(input=bam)
    await assert_md5(out, expected="<md5 from main.nf.test.snap>", label="samtools sort")
```

Then export it from the module's `tests` list and (if it's in a new file)
import that list in [`tests/run_all.py`](tests/run_all.py):

```python
# tests/modules/test_samtools.py
tests = [test_samtools_sort, ...]

# tests/run_all.py
from tests.modules import test_samtools
ALL_TESTS = [*test_bedtools.tests, *test_samtools.tests, ...]
```

### 4. Run it

Tests run against an actual Flyte cluster via the native `flyte run` CLI.
There's no custom entry point or `--mode` flag — the target cluster,
project, domain, and auth all come from your Flyte config (`--config`
flag, `FLYTE_CONFIG`, or `~/.flyte/config.yaml`). Run from the project
root so the `tests` package is importable.

**Devbox** (the fast inner loop). First, in another terminal:

```bash
flyte start devbox
```

That spins up a local Union cluster at `localhost:30080` and points your
default config at it. Then:

```bash
flyte run tests/run_all.py run_all
```

Each test gets its own pod on the devbox, the `fixture` task is
cached so each dataset downloads once, and the wrapped tool tasks reuse
their biocontainer images across actions. Because `flyte run` reads the
config (a `localhost` endpoint) before building the image, the test image
auto-targets the devbox's local registry — no extra flags.

**Remote** (full CI run, scale horizontally) — point at another config:

```bash
flyte --config ~/.flyte/prod.yaml run tests/run_all.py run_all
```

Same test code, same Flyte caching semantics — just a different config.
The default push registry for a non-devbox cluster is `ghcr.io/flyteorg`
(which only Flyte CI can push to); if your cluster's registry differs,
set `FLYTE_IMAGE_REGISTRY` to it. To run without any cluster, use
`flyte run --local tests/run_all.py run_all`.

To iterate on a subset, write a small sibling driver that imports only
the tests you care about:

```python
# tests/modules/run_bedtools.py
from tests.framework import env, gather_tests
from tests.modules.test_bedtools import tests

@env.task
async def run() -> str:
    return await gather_tests(tests)
```

Invoke it the same way: `flyte run tests/modules/run_bedtools.py run`.

If an md5 assert fails, the summary lists which test failed and the
expected-vs-actual hex. The wrapped task's stdout/stderr show up in the
Flyte UI for that action's pod.
