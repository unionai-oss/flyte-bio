# Contributing

flyte-bio's modules and pipelines are ports of
[nf-core](https://nf-co.re/) modules and pipelines: same containers, same
command lines, same test data. Following upstream that closely keeps the
behaviour familiar. It also gives us a regression oracle: each upstream
module ships test cases with snapshot md5s we can assert against.

## Porting a module

### 1. Find the upstream version

Port the module **at the version the target pipeline pins**, not upstream's
latest. The pipeline's `modules.json` gives each module's `git_sha`. Take the
container and the snapshot from that same SHA, so the md5s match the
container:

```bash
SHA=$(curl -s https://raw.githubusercontent.com/nf-core/sarek/3.10.0/modules.json \
  | jq -r '.repos["https://github.com/nf-core/modules.git"].modules["nf-core"]["bwa/index"].git_sha')
for f in main.nf tests/main.nf.test tests/main.nf.test.snap tests/nextflow.config; do
  curl -s https://raw.githubusercontent.com/nf-core/modules/$SHA/modules/nf-core/bwa/index/$f
done
```

What to take from each:

- **`main.nf`:** the container and the command line.
  - Use the Docker branch of the `container` line, not the Singularity URL.
    That's usually `community.wave.seqera.io/...`. A bare
    `biocontainers/...` resolves to Docker Hub, which the cluster can't pull
    from, so rewrite it to its `quay.io/biocontainers/...` mirror.
  - The `process_*` label hints at resources.
- **The pipeline's `conf/modules/*.config`:** the `ext.args` and `ext.prefix`
  the pipeline uses. They often matter as much as the module itself.
- **`tests/main.nf.test`, `tests/nextflow.config`:** the test cases, their
  inputs (paths under
  [nf-core/test-datasets@modules](https://github.com/nf-core/test-datasets/tree/modules)),
  and any `ext.args` the tests set.
- **`tests/main.nf.test.snap`:** the expected outputs. Never use a `- stub`
  case's md5s: stub runs only `touch` their outputs.

### 2. Write the shell task

One file per tool family under `src/flyte_bio/modules/`. Each file has a
module docstring that describes the tools and their output names. Follow an
existing module of the same shape: `bwa.py` for a simple tool,
`gatk4.py` for several tools with multiple images, `umitools.py` for tools
with an output directory and a wrapper.

```python
"""bwa — BWA short-read alignment. ..."""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BWA_IMAGE = "community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0"
INDEX_RESOURCES = flyte.Resources(cpu=1, memory="24Gi")

bwa_index_cmd = shell.create(
    name="bwa_index",
    image=BWA_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=r"""
        FA=({inputs.fasta})
        NAME=$(basename "${FA[0]}")
        ARGS={inputs.args}
        bwa index $ARGS -p {outputs.index}/"${NAME%.*}" "${FA[0]}"
    """,
)

env = flyte.TaskEnvironment.from_task("bwa", bwa_index_cmd.as_task())


async def bwa_index(fasta: File, args: str = "") -> Dir:
    """BWA index of ``fasta``: a Dir of ``<fasta name>.{amb,ann,bwt,pac,sa}``."""
    return await bwa_index_cmd(fasta=fasta, args=args)
```

Rules that keep tasks correct:

- **Inputs:**
  - **Expand every file input into an array first:** `FA=({inputs.fasta})`,
    then use `"${FA[0]}"`. A `File` input renders as a glob over its
    staging directory, so the file keeps its original name.
  - **Optional inputs** (`File | None`, `list[File]`) expand under
    `shopt -s nullglob`, so an unset input is an empty array.
  - **`args` is expanded unquoted** (`ARGS={inputs.args}`, then `$ARGS`), so
    it can carry several options. If options can contain regex characters
    (`*`, `?`), run `set -f` after expanding the file globs.
- **Outputs:**
  - **When there are several, write them under their upstream names into one
    output `Dir`.** The wrapper then picks files out with
    `results.get_file(name)`. A bare `File` output is named after its output
    key and has no extension, which breaks tools that sniff extensions.
  - **Tools that look for companion files** (an index next to its BAM, a
    `.fai` and `.dict` next to the FASTA, a `.tbi` next to its VCF) need them
    symlinked side by side in a working directory: `W=$(mktemp -d)`. Inputs
    are read-only, and some tools (htslib, samtools) write a `.fai` next to
    the FASTA.
  - **A script whose last command is `[ test ] && cmd`** exits non-zero when
    the test is false; end it with `true`.
- **Environments:** one image per `TaskEnvironment`. A family whose tools use
  different images has one environment per image and an aggregate `env` (see
  `gatk4.py`). Add the module's `env` to `flyte_bio.modules.env` in
  `src/flyte_bio/modules/__init__.py`, and its tests to `tests/run_all.py`.
- **Helper scripts** that upstream ships in `bin/` or `templates/` (Python,
  R, Perl) are vendored unchanged under `src/flyte_bio/scripts/`:
  - Undo Groovy's escaping (`\\` → `\`, `\$` → `$`).
  - Turn templated values into CLI arguments.
  - Register the file in `ship_in_bundle(...)`.
  - The wrapper passes the script to a stock interpreter image as a `File`
    input (see `tx2gene.py`). Don't rewrite a helper in another language.
- **Naming:** name things by what they do, not by the upstream project.
  Provenance goes in docs and functional URLs.

### 3. Test it

Port the upstream module's test cases into `tests/modules/test_<module>.py`.
Tests are Flyte tasks (below). Which assertion to use depends on whether the
upstream md5 can be reproduced:

| Output | Assertion |
|---|---|
| Pure data: tables, indexes, counts, BED/VCF bodies without dates | `assert_md5(file, md5)` |
| `.gz` outputs (nf-test hashes them decompressed) | `assert_gunzipped_md5` |
| BAM/CRAM where upstream snapshots `getReadsMD5` / `md5Reads` | `assert_reads_md5` (the read sequences' md5) |
| Files whose header records a command line, path or date: BAM headers, samtools stats, most VCFs, MarkDuplicates metrics | content checks: the expected header, record counts, invariants (e.g. BQSR doesn't change sequences) |

Thread counts can change outputs: for example, featureCounts' header records
`-T`. When an md5 doesn't match, check upstream's test `ext.args` and thread
counts before suspecting the tool.

## Porting a pipeline

A pipeline is a module under `src/flyte_bio/pipelines/` made of plain async
functions, not tasks. They run in the caller's task and call module tasks.
To match the existing pipelines:

- **One public function per upstream subworkflow** (`prepare_reference`,
  `map_reads`, …), plus a main function that takes the samplesheet `File`,
  like `--input`, and the upstream parameters as keyword arguments.
- **Read the samplesheet with `flyte_bio.samplesheet.read_samplesheet`** and
  index the upstream columns directly. Validate everything before starting
  work.
- **Return a result dataclass** of `File`s and `Dir`s. Write a `*_layout()`
  mapping it onto upstream's `--outdir` folders, and publish with
  `flyte_bio.publish.publish` when `outdir` or `publish_results` is set.
- **Keep upstream's file names** (prefixes from the pipeline's
  `conf/modules`): MultiQC and users both depend on them.
- **Pipeline tests** (`tests/pipelines/`) run upstream's test profiles end to
  end. List every divergence from upstream in the pipeline's module docstring
  and its page under `docs/pipelines/`.

## Running the tests

Tests are Flyte tasks, so they run on a cluster like any workflow, with one
pod per test. `tests/run_all.py` gathers them all. The cluster, project and
auth come from your Flyte config:

```bash
flyte run tests/run_all.py run_all                          # default config
flyte --config ~/.flyte/other.yaml run tests/run_all.py run_all
flyte run --local tests/run_all.py run_all                  # no cluster
```

Run from the repository root. The test image installs this checkout
(`with_uv_project`), so the code under test is what's on disk.

To iterate on a subset, write a small driver next to the tests:

```python
# tests/modules/run_bwa.py
from tests.framework import env, gather_tests
from tests.modules.test_bwa import tests


@env.task
async def run() -> str:
    return await gather_tests(tests)
```

```bash
flyte run tests/modules/run_bwa.py run
```

Helpers in [`tests/framework.py`](../tests/framework.py):

- **Fixtures:**
  - `fixture(path)` fetches a test-data file, by a path under the
    test-datasets `modules` branch or by absolute URL. It's cached, so each
    file downloads once.
  - `fixture_dir(url, member)` does the same for a `.tar.gz` fixture,
    extracted to a `Dir`.
  - `stage_samplesheet(url, file_columns)` copies an upstream samplesheet,
    swapping its file URLs for cached fixtures.
- **Assertions:** `assert_md5`, `assert_gunzipped_md5`, `assert_reads_md5`
  (with `reads_md5`), `assert_nonempty`, `assert_dir_nonempty`.
- **Running:** `gather_tests(tests)` runs tests in parallel, prints a summary
  and fails if any test failed.

### Debugging a failing tool

The shell wrapper sends a tool's stdout and stderr to files that aren't
uploaded when the task fails. To see them, rerun the task with the same
inputs and its script wrapped as `{ ...; } 2>&1 | tee /proc/1/fd/2`. That
sends the output to the container log, which `flyte get logs <run>` shows.

### Cached results after editing a script

Shell tasks are cached on their inputs, and today every shell task gets the
same cache version, so editing a module's script doesn't invalidate its cache
(fix in flyteorg/flyte-sdk#1657). Until that lands, after changing a script,
rename the task or set `cache="disable"` while you iterate.
