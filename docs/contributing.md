# Contributing

flyte-bio's tools are Flyte shell tasks that run each tool's command line in
a pinned container, and its workflows are plain async Python that composes
them. Every tool is tested against known-good outputs.

## Writing a tool

One file per tool family under `src/flyte_bio/modules/`. Each file has a
module docstring that describes the tools and their output names. Start from
an existing module of the same shape: `bwa.py` for a simple tool,
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
  - **When there are several, write them under their conventional names into
    one output `Dir`.** The wrapper then picks files out with
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
- **Helper scripts** (Python, R, Perl) that a container doesn't include are
  vendored unchanged under `src/flyte_bio/scripts/`:
  - Register the file in `ship_in_bundle(...)`.
  - The wrapper passes the script to a stock interpreter image as a `File`
    input (see `tx2gene.py`). Don't rewrite a helper in another language.
- **Threads and memory** come from the container's limits at run time:
  start the script with `{TASK_LIMITS}` and use `$TASK_CPUS` and
  `$TASK_MEMORY_MB` (see [`limits.py`](../src/flyte_bio/limits.py)), so a
  resource override reaches the tool.
- **Exit codes:** end commands that can exit 255 with `|| exit 1`; Flyte
  retries exit code 255 as a system error.
- **Naming:** name things by what they do. Provenance goes in docs and
  functional URLs.
- **Licences:** tools with restrictive licences (EULAs, non-commercial,
  vendor-only) go in `licensed_env` in `src/flyte_bio/modules/__init__.py`,
  with the reason in the comment above it, and the licence named in their
  row of the [modules table](modules.md#available-modules).
- **Docs:** add a row to the [modules table](modules.md#available-modules).

## Testing a tool

Tests live in `tests/modules/test_<module>.py` as Flyte tasks (below), with a
module-level `tests` list that `tests/run_all.py` collects. Each test runs
the tool on small fixture data and checks its outputs:

| Output | Assertion |
|---|---|
| Pure data: tables, indexes, counts, BED/VCF bodies without dates | `assert_md5(file, md5)` |
| `.gz` outputs | `assert_gunzipped_md5` (md5 of the decompressed bytes) |
| BAM/CRAM | `assert_reads_md5` (the read sequences' md5) |
| Files whose header records a command line, path or date: BAM headers, samtools stats, most VCFs | content checks: the expected header, record counts, invariants |
| Floating-point values whose trailing digits vary by CPU | value checks within a tolerance (see `tests/modules/test_dupradar.py`) |

Thread counts can change outputs: for example, featureCounts' header records
`-T`. Pin such a test's CPUs with `tool.override(resources=...)` and a cache
salt. Test task names must be unique across the suite: prefix them
`test_<module>_`. A test that can't run yet (licensed data, missing
fixtures) raises `NotTested(reason)` and is reported as skipped.

Bringing a tool over from Nextflow? [Migrating from
Nextflow](migrating-from-nextflow.md) covers finding the upstream version and
reusing its test snapshots.

## Writing a workflow

A workflow (`src/flyte_bio/pipelines/`) is a module of plain async
functions, not tasks: they run in the caller's task and call tool tasks.

- **One public function per stage** (`prepare_reference`, `map_reads`, …),
  plus a main function that takes the samplesheet `File` and the options as
  keyword arguments.
- **Read the samplesheet with `flyte_bio.samplesheet.read_samplesheet`** and
  validate everything before starting work.
- **Return a result dataclass** of `File`s and `Dir`s. Write a `*_layout()`
  mapping it onto an output folder tree, and publish with
  `flyte_bio.publish.publish` when `outdir` or `publish_results` is set.
- **Add a Hydra config** under `src/flyte_bio/conf/<workflow>/`.
- **Tests** (`tests/pipelines/`) run the workflow end to end on a small test
  profile.

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
  - `fixture(path)` fetches a test-data file, by a path under the shared
    test-datasets repository or by absolute URL. It's cached, so each
    file downloads once.
  - `fixture_dir(url, member)` does the same for a `.tar.gz` fixture,
    extracted to a `Dir`.
  - `stage_samplesheet(url, file_columns)` copies a samplesheet,
    swapping its file URLs for cached fixtures.
- **Assertions:** `assert_md5`, `assert_gunzipped_md5`, `assert_reads_md5`
  (with `reads_md5`), `assert_nonempty`, `assert_dir_nonempty`.
- **Running:** `gather_tests(tests)` runs tests in parallel, prints a summary
  (passed, failed, not tested) and fails if any test failed.
- **Skipping:** `NotTested(reason)` marks a test that can't run yet.

### Debugging a failing tool

The shell wrapper sends a tool's stdout and stderr to files that aren't
uploaded when the task fails. To see them, rerun the task with the same
inputs and its script wrapped as `{ ...; } 2>&1 | tee /proc/1/fd/2`. That
sends the output to the container log, which `flyte get logs <run>` shows.

### Cached results after editing a script

Shell tasks are cached on their inputs and versioned by their image, command
and arguments (flyte 2.12 and later), so editing a module's script reruns it.
