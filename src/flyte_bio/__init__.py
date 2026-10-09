"""flyte-bio — curated bioinformatics building blocks for Flyte.

Two sub-packages:

- :mod:`flyte_bio.modules` — typed shell-task wrappers around individual
  bio CLIs (bedtools, samtools, ...). Each module exposes its tasks and a
  module-level ``env`` (a :class:`flyte.TaskEnvironment`).
- :mod:`flyte_bio.pipelines` — higher-level pipelines composed from the
  module tasks (rnaseq, variant calling, ...).

For convenience this package re-exports :data:`flyte_bio.modules.env`, an
aggregate environment that depends on every module env. Pipelines can use
it as a single ``depends_on`` entry to pull in the whole module suite::

    import flyte
    from flyte_bio import env as bio_env

    env = flyte.TaskEnvironment(
        name="genomics_pipeline",
        depends_on=[bio_env],
    )
"""



from .modules import env

__all__ = ["env"]
