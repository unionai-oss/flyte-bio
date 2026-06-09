"""flyte_bio.pipelines — higher-level pipelines composed from module tasks.

Empty for now. Pipelines like ``rnaseq`` and variant-calling workflows will
live as submodules here (e.g. ``flyte_bio.pipelines.rnaseq``), each
exposing its own ``TaskEnvironment`` that depends on whichever module envs
its tasks need.
"""
