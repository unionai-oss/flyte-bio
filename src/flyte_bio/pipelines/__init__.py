"""flyte_bio.pipelines — higher-level pipelines composed from module tasks.

- :mod:`flyte_bio.pipelines.rnaseq` — RNA-seq STAR alignment + salmon
  quantification (the ``star_salmon`` path).

Pipelines are plain async functions that run inside the caller's task and
fan out to module tasks, so the caller's ``TaskEnvironment`` must
``depends_on`` :data:`flyte_bio.modules.env`.
"""
