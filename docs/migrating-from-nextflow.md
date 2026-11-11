# Migrating a Nextflow module or pipeline

Most of flyte-bio's tools started life as nf-core modules, and its workflows
as nf-core pipelines. This guide covers bringing another one across: where to
find the upstream module at the right version, how its parts map onto a Flyte
shell task, and how to reuse its test snapshots as a regression oracle.

For the general rules on writing a tool and its tests, which apply whatever
the source, see [contributing](contributing.md).

## Concept map

| Nextflow / nf-core | flyte-bio |
|---|---|
| `process` in `modules/nf-core/<tool>/<subtool>/main.nf` | a shell task (`flyte.extras.shell.create`) in `src/flyte_bio/modules/<tool>.py` |
| `container` | the shell task's `image` (Docker branch, not Singularity) |
| `label 'process_*'` | the task's `flyte.Resources` |
| `task.cpus`, `task.memory` | `$TASK_CPUS`, `$TASK_MEMORY_MB` from `{TASK_LIMITS}` ([`limits.py`](../src/flyte_bio/limits.py)) |
| `ext.args` | the `args` input |
| `ext.prefix` | the `prefix` input |
| `bin/`, `templates/` scripts | vendored unchanged under `src/flyte_bio/scripts/` |
| `tests/main.nf.test` + `.snap` | `tests/modules/test_<tool>.py` with the snapshot md5s |
| subworkflow | a public async function in `src/flyte_bio/pipelines/<pipeline>.py` |
| `--input` samplesheet | the pipeline function's samplesheet `File` |
| `--outdir` | `outdir=` / `publish_results=True` and a `*_layout()` function |
| `-params-file` | an OmegaConf config ([`flyte_bio.config`](../src/flyte_bio/config.py)) |
| `baseName` / `simpleName` | a file's name without its last / any extension |

## 1. Find the upstream version

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
    `biocontainers/...` resolves to Docker Hub, which clusters often can't
    pull from, so rewrite it to its `quay.io/biocontainers/...` mirror.
  - The `process_*` label hints at resources.
- **The pipeline's `conf/modules/*.config`:** the `ext.args` and `ext.prefix`
  the pipeline uses. They often matter as much as the module itself.
- **`tests/main.nf.test`, `tests/nextflow.config`:** the test cases, their
  inputs (paths under
  [nf-core/test-datasets@modules](https://github.com/nf-core/test-datasets/tree/modules),
  which `tests.framework.fixture` resolves), and any `ext.args` the tests set.
- **`tests/main.nf.test.snap`:** the expected outputs. Never use a `- stub`
  case's md5s: stub runs only `touch` their outputs.

If a module at the pinned SHA is deprecated or a stub, port its last working
version and record the commit in the module's docstring.

## 2. Translate the process

Write the shell task as described in
[contributing](contributing.md#writing-a-tool). Nextflow-specific points:

- **Undo Groovy's escaping** in scripts and templates: `\\` → `\`,
  `\$` → `$`. Templated values (`${meta.id}`, `${task.cpus}`) become CLI
  arguments or shell task inputs.
- **Outputs keep upstream's names** (`<prefix>.md.cram`,
  `<prefix>.md.cram.metrics`), written into one output `Dir` when there are
  several. MultiQC and users depend on them.
- **Staging:** Nextflow stages inputs by name into the work directory; a
  shell task stages each input under its own name in a read-only directory.
  Tools that need companion files side by side, or write next to their
  inputs, need a scratch working directory (`W=$(mktemp -d)`) with symlinks.
- **Name things by what they do**, not by the upstream project. Provenance
  belongs in docs and functional URLs.

## 3. Port the tests

Port every non-stub upstream test case into `tests/modules/test_<tool>.py`,
with a `# upstream case: <tool/subtool> "<test name>"` comment. Choose the
assertion by whether upstream's md5 can be reproduced:

| Output | Assertion |
|---|---|
| Pure data: tables, indexes, counts, BED/VCF bodies without dates | `assert_md5(file, md5)` |
| `.gz` outputs (nf-test hashes them decompressed) | `assert_gunzipped_md5` |
| BAM/CRAM where upstream snapshots `getReadsMD5` / `md5Reads` | `assert_reads_md5` (the read sequences' md5) |
| Files whose header records a command line, path or date: BAM headers, samtools stats, most VCFs, MarkDuplicates metrics | content checks: the expected header, record counts, invariants |
| Floating-point values whose trailing digits vary by CPU | value checks within a tolerance (see `tests/modules/test_dupradar.py`) |

Thread counts can change outputs: for example, featureCounts' header records
`-T`. When an md5 doesn't match, check upstream's test `ext.args` and thread
counts before suspecting the tool, and pin the test's CPUs if needed.

## Migrating a pipeline

A pipeline becomes a module under `src/flyte_bio/pipelines/` made of plain
async functions, not tasks. They run in the caller's task and call tool tasks.

- **One public function per upstream subworkflow** (`prepare_reference`,
  `map_reads`, …), plus a main function that takes the samplesheet `File` and
  the upstream parameters as keyword arguments.
- **Port the pipeline's local modules** too: they aren't in the shared
  modules repository.
- **Read the samplesheet with `flyte_bio.samplesheet.read_samplesheet`** and
  index the upstream columns directly. Validate everything before starting
  work.
- **Return a result dataclass** of `File`s and `Dir`s. Write a `*_layout()`
  mapping it onto upstream's `--outdir` folders, and publish with
  `flyte_bio.publish.publish` when `outdir` or `publish_results` is set.
- **Keep upstream's file names** (prefixes from the pipeline's
  `conf/modules`).
- **Add a Hydra config** under `src/flyte_bio/conf/<pipeline>/`, with
  `profile/test.yaml` pointing at upstream's test profile data.
- **Pipeline tests** (`tests/pipelines/`) run upstream's test profiles end to
  end. List every divergence from upstream in the pipeline's module docstring
  and its page under `docs/pipelines/`. Report-level differences are expected:
  MultiQC sections that only Nextflow produces (run parameters, software
  versions, methods) are left out.

## Coverage

[`porting-plan.md`](porting-plan.md) explains how upstream modules were
ranked for porting, and [`module-popularity.csv`](module-popularity.csv)
tracks every module used by a released nf-core pipeline, with whether it's
ported.
