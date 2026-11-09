"""gnusort — sorting text files (BED, GTF, CSV, ...) with GNU coreutils ``sort``.

Exposes :func:`gnu_sort`, which runs ``sort <args> <input>`` and writes the
result to ``<prefix>.<suffix>`` (upstream's prefix defaults to the sample
id). ``args`` are sort's options, e.g. ``-k1,1 -k2,2n`` for a BED or
``-k1,1 -k4,4n`` for a GTF; a field separator goes in bare (``-t ;``), as
``args`` is split on whitespace and not otherwise parsed by a shell. As
upstream, the container's locale is left as it is.
"""

import os

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

COREUTILS_IMAGE = "quay.io/biocontainers/coreutils:9.5"
# Upstream labels it process_low; sort's buffer is left at its default.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

gnu_sort_cmd = shell.create(
    name="gnu_sort",
    image=COREUTILS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"input": File, "suffix": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.input})
        S={inputs.suffix}
        P={inputs.prefix}
        ARGS={inputs.args}
        set -f
        sort $ARGS "${IN[0]}" > {outputs.results}/"$P.$S"
    """,
)

env = flyte.TaskEnvironment.from_task("gnusort", gnu_sort_cmd.as_task())


async def gnu_sort(input: File, suffix: str, prefix: str = "sample", args: str = "") -> File:
    """``input`` sorted by ``sort <args>``, as ``<prefix>.<suffix>``."""
    output = f"{prefix}.{suffix}"
    if (input.name or os.path.basename(input.path)) == output:
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    results = await gnu_sort_cmd(input=input, suffix=suffix, prefix=prefix, args=args)
    sorted_file = await results.get_file(output)
    if sorted_file is None:
        raise FileNotFoundError(f"sort wrote no {output}")
    return sorted_file
