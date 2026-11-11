"""flyte-bio — the Flyte 2 open-source library for biology.

Two sub-packages:

- :mod:`flyte_bio.modules` — scientific tools (bedtools, samtools, GATK, ...)
  as typed, containerized Flyte tasks. Each module exposes its tasks and a
  module-level ``env`` (a :class:`flyte.TaskEnvironment`).
- :mod:`flyte_bio.pipelines` — end-to-end workflows composed from those
  tools (rnaseq, scrnaseq, variant calling).

For convenience this package re-exports :data:`flyte_bio.modules.env`, an
aggregate environment that depends on every module env. A workflow can use
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
