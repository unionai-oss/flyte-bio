# Known issues

Temporary workarounds for bugs outside flyte-bio. Each says what to undo once
the fix is available.

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
