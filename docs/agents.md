# flyte-bio for agents

flyte-bio is the Flyte 2 open-source library for biology: about 900
scientific tools as typed, containerized Flyte tasks (`flyte_bio.modules`),
and end-to-end workflows built from them (`flyte_bio.pipelines`).

## Using the library

- **Find a tool:** search the catalog in `docs/modules.md` (one row per
  family: what it does, function names). Then read the module's docstring in
  `src/flyte_bio/modules/<family>.py` for inputs, outputs and file names.
- **Call a tool** from your own `@env.task` with `await`, e.g.
  `await samtools_sort(bam=bam)`. Data is `flyte.io.File` / `Dir`; options
  are plain values; extra tool flags go in `args` (whitespace-split).
  Multi-output tools return a small dataclass.
- **Your environment** needs flyte-bio in its image and
  `depends_on=[flyte_bio.env]`. Licensed tools (vendor EULAs, non-commercial
  terms) are in `flyte_bio.modules.licensed_env`, which `flyte_bio.env`
  leaves out: add it only after the user has accepted the terms.
- **Resources:** `await tool.override(resources=flyte.Resources(cpu=8, memory="32Gi"))(...)`.
  The tool's thread and memory flags follow the resources it gets.
- **Workflows:** `rnaseq`, `scrnaseq`, `variant_calling` take a samplesheet
  `File`; see `docs/pipelines/`. Runnable scripts are in `examples/`.

## Changing the library

- Read `docs/contributing.md` first. Bringing a tool over from Nextflow:
  `docs/migrating-from-nextflow.md`.
- Shell scripts start with `{TASK_LIMITS}` and take threads and memory from
  `$TASK_CPUS` / `$TASK_MEMORY_MB`. End commands that can exit 255 with
  `|| exit 1`.
- Register a new module's `env` in `src/flyte_bio/modules/__init__.py`, its
  tests in `tests/run_all.py`, and a row in `docs/modules.md`.
- Tests are Flyte tasks, run on a cluster:
  `flyte run tests/run_all.py run_all`, or a small driver over one module's
  `tests` list. `NotTested(reason)` marks a test that can't run yet.
- Style: Python 3.12+, ruff (line length 120); no
  `from __future__ import annotations`; no leading-underscore names; name
  things by what they do.
