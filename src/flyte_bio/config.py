"""Configure a pipeline from an OmegaConf config, like Nextflow's ``-params-file``.

A pipeline's config is flat: one key per keyword argument of the pipeline
function, with the same name and default. Files and directories are given as
URIs (``s3://…``, ``gs://…``, ``https://…``)::

    protocol: 10XV3
    fasta: s3://my-bucket/ref/genome.fa
    gtf: s3://my-bucket/ref/genes.gtf
    skip_cellbender: true

Run it from your own task, which takes the config as a ``DictConfig``
(this needs ``flyte-bio[config]`` in the task's image)::

    @env.task
    async def run(cfg: DictConfig, samplesheet: File | None = None) -> Dir:
        result = await run_pipeline(scrnaseq, cfg, samplesheet=samplesheet)
        return result.outdir

Keyword arguments passed to :func:`run_pipeline` take precedence over the
config. They are how a local file reaches the pipeline: ``flyte run`` and
``flyte hydra run`` upload a task's ``File`` arguments, but a local path in
a config is just a string.
"""

import difflib
import inspect
import types
import typing
from collections.abc import Awaitable, Callable, Mapping
from typing import Any, TypeVar

import flyte
from flyte.io import Dir, File
from omegaconf import MISSING, DictConfig, OmegaConf

from flyte_bio.conf import CONFIG_DIR

R = TypeVar("R")

__all__ = ["CONFIG_DIR", "bind_config", "pipeline_defaults", "run_pipeline"]

# Keys that Hydra and its Flyte launcher put in a composed config.
RESERVED_KEYS = frozenset({"hydra", "task_env"})


def pipeline_defaults(pipeline: Callable[..., Any]) -> DictConfig:
    """The config of ``pipeline``'s defaults; its required options are ``???`` (missing)."""
    values = {}
    for name, param in inspect.signature(pipeline).parameters.items():
        if param.default is inspect.Parameter.empty:
            values[name] = MISSING
        elif isinstance(param.default, tuple):
            values[name] = list(param.default)
        else:
            values[name] = param.default
    return OmegaConf.create(values)


def bind_config(pipeline: Callable[..., Any], cfg: Mapping[str, Any] | DictConfig, **given: Any) -> dict[str, Any]:
    """Check ``cfg`` against ``pipeline``'s signature and convert it into its keyword arguments.

    ``given`` arguments are used as they are, and take precedence over ``cfg``.
    Raises ValueError for an unknown option, a missing required one, or a value
    of the wrong type.
    """
    params = inspect.signature(pipeline).parameters
    hints = typing.get_type_hints(pipeline)
    values = OmegaConf.to_container(cfg, resolve=True) if isinstance(cfg, DictConfig) else dict(cfg)
    assert isinstance(values, dict)

    unknown = [k for k in values if k not in params and k not in RESERVED_KEYS]
    unknown += [k for k in given if k not in params]
    if unknown:
        hints_text = []
        for key in unknown:
            close = difflib.get_close_matches(key, params, n=1)
            hints_text.append(f"{key!r}" + (f" (did you mean {close[0]!r}?)" if close else ""))
        raise ValueError(f"unknown {pipeline.__name__} options: {', '.join(hints_text)}")

    kwargs = {}
    errors = []
    for name, param in params.items():
        if given.get(name) is not None:
            kwargs[name] = given[name]
            continue
        value = values.get(name, MISSING)
        if value is MISSING or value == "???":
            if param.default is inspect.Parameter.empty:
                errors.append(f"{name}: required")
            continue
        try:
            kwargs[name] = convert(value, hints[name])
        except (TypeError, ValueError) as e:
            errors.append(f"{name}: {e}")
    if errors:
        raise ValueError(f"invalid {pipeline.__name__} config:\n  " + "\n  ".join(errors))
    return kwargs


async def run_pipeline(pipeline: Callable[..., Awaitable[R]], cfg: Mapping[str, Any] | DictConfig, **given: Any) -> R:
    """Run ``pipeline`` with the options in ``cfg``, and ``given`` arguments over them."""
    return await pipeline(**bind_config(pipeline, cfg, **given))


def convert(value: Any, annotation: Any) -> Any:
    """``value`` (a plain config value) as an argument of type ``annotation``."""
    origin = typing.get_origin(annotation)
    args = typing.get_args(annotation)

    if origin in (typing.Union, types.UnionType):
        if value is None:
            if type(None) in args:
                return None
            raise TypeError("can't be null")
        (inner,) = [a for a in args if a is not type(None)]
        return convert(value, inner)
    if value is None:
        raise TypeError("can't be null")

    if annotation in (File, Dir):
        return location(value, annotation)
    if origin is list:
        (item,) = args
        if not isinstance(value, list):
            raise TypeError(f"expected a list, got {value!r}")
        return [convert(v, item) for v in value]
    if origin is tuple:
        if not isinstance(value, list):
            raise TypeError(f"expected a list, got {value!r}")
        return tuple(convert(v, args[0]) for v in value)
    if origin is dict:
        key_type, item = args
        if not isinstance(value, dict):
            raise TypeError(f"expected a mapping, got {value!r}")
        return {convert(k, key_type): convert(v, item) for k, v in value.items()}
    if annotation is float and isinstance(value, int) and not isinstance(value, bool):
        return float(value)
    if annotation in (bool, int, float, str):
        # bool is a subclass of int, so check it separately.
        if isinstance(value, annotation) and (annotation is bool or not isinstance(value, bool)):
            return value
        raise TypeError(f"expected {annotation.__name__}, got {value!r}")
    raise TypeError(f"unsupported option type {annotation}")


def location(value: Any, kind: type[File] | type[Dir]) -> File | Dir:
    """The ``File`` or ``Dir`` at ``value``, a URI (or, outside a cluster, a local path)."""
    if not isinstance(value, str):
        raise TypeError(f"expected a URI, got {value!r}")
    ctx = flyte.ctx()
    if "://" not in value and ctx is not None and ctx.is_in_cluster():
        raise ValueError(
            f"{value!r} is a local path, which the cluster can't read: give a URI, "
            "or pass the file as a task argument so flyte run uploads it"
        )
    return kind.from_existing_remote(value)
