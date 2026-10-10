# Known issues

Temporary workarounds for bugs outside flyte-bio. Each says what to undo once
the fix is available.

## MultiQC plot export is off

**Affects:** both pipelines' MultiQC reports
(`scripts/multiqc_config.yml` and `scripts/multiqc_config_variant_calling.yml`).

Both vendored MultiQC configs set `export_plots: false` (upstream: `true`).
The MultiQC output keeps the HTML report and data tables, but not the static
PNG/SVG/PDF plot exports.

**Why:** with plot exports on, the output is a few hundred files, and the
Flyte copilot uploader sidecar runs out of memory uploading them. The task
fails with "Pod reported success despite being OOMKilled". Two things combine:
- since stow v0.5.0 (flyteorg/flyte#8115), every uploaded file is read into
  memory;
- copilot uploads all of a directory's files at once, and the sidecar's
  memory limit is 128Mi.

**Undo** by setting both configs back to `export_plots: true` once either
fix is merged **and deployed**. On Union clusters, deployed means the
`flyte2` submodule pin in the `cloud` repo includes the fix and the cluster
has been redeployed.

- [flyteorg/flyte#8149](https://github.com/flyteorg/flyte/pull/8149): copilot
  uploads at most 8 files at once.
- [flyteorg/stow#32](https://github.com/flyteorg/stow/pull/32): stow streams
  file bodies instead of buffering them. This one also needs a stow release
  and a `go.mod` bump in flyte.

## Editing a module's script doesn't invalidate its cache

Every shell task currently gets the same cache version, so after you change a
module's script, cached results from the old script can still be used. Rename
the task or disable its cache while iterating.

**Fix:** [flyteorg/flyte-sdk#1657](https://github.com/flyteorg/flyte-sdk/pull/1657),
which versions container tasks by image, command and arguments.

## Shell tasks can't override resources per call

`task.override(resources=...)` fails on shell tasks with
`unexpected keyword argument 'short_name'`, a flyte-sdk bug. Resources are
set per module (constants near the top of each module file).

## `https://` files can't be read inside a task

Reading an `https://` `File` from inside a task, e.g. a samplesheet, fails
with an aiohttp event-loop error. Shell tasks fetch `https://` inputs fine,
so this only affects files the pipeline code reads itself:
- **samplesheets:** pass them as `s3://` / `gs://` or as a local file to
  `flyte run`;
- **FASTQ headers**, which variant_calling reads with the standard library
  instead.
