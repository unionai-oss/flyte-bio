"""wget — download a file over HTTP(S) or FTP with GNU Wget.

Exposes :func:`wget`, which downloads ``url`` inside the task and returns it as
``<prefix>.<suffix>``, as upstream (``wget -O -`` to that file). ``args`` are
extra ``wget`` options.

The task is cached on the URL, the name and ``args``, so a URL downloads once.
Raise its disk with ``wget_cmd.override(resources=...)`` for large
files.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

WGET_IMAGE = "community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi", disk="10Gi")

wget_cmd = shell.create(
    name="wget",
    image=WGET_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"url": str, "prefix": str, "suffix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        URL={inputs.url}
        NAME={inputs.prefix}.{inputs.suffix}
        ARGS={inputs.args}
        OUT={outputs.results}

        wget \
            -O - \
            $ARGS \
            "$URL" \
            > "$OUT/$NAME" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("wget", wget_cmd.as_task())


async def wget(url: str, suffix: str, prefix: str = "sample", args: str = "") -> File:
    """Download ``url`` as ``<prefix>.<suffix>``."""
    results = await wget_cmd(url=url, prefix=prefix, suffix=suffix, args=args)
    name = f"{prefix}.{suffix}"
    downloaded = await results.get_file(name)
    if downloaded is None:
        raise FileNotFoundError(f"wget downloaded no {name}")
    return downloaded
