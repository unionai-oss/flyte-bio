"""aria2 — download a file over HTTP(S), FTP, BitTorrent or Metalink with aria2c.

Exposes :func:`aria2`, which downloads ``source_url`` inside the task and
returns the downloaded file under its name in the URL (the URL's last path
component), as upstream. ``args`` are extra ``aria2c`` options (for example
``-x 4`` for more connections per server). As upstream, certificates aren't
checked (``--check-certificate=false``).

The task is cached on the URL and ``args``, so a URL downloads once. Raise its
disk with ``aria2_cmd.override(resources=...)`` for large files.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

ARIA2_IMAGE = "community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi", disk="10Gi")

aria2_cmd = shell.create(
    name="aria2",
    image=ARIA2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"source_url": str, "file_name": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        URL={inputs.source_url}
        NAME={inputs.file_name}
        ARGS={inputs.args}
        OUT={outputs.results}

        aria2c \
            --check-certificate=false \
            --dir "$OUT" \
            --out "$NAME" \
            $ARGS \
            "$URL" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("aria2", aria2_cmd.as_task())


async def aria2(source_url: str, args: str = "") -> File:
    """Download ``source_url``: the file, named after the URL's last path component."""
    name = source_url.rsplit("/", maxsplit=1)[-1]
    results = await aria2_cmd(source_url=source_url, file_name=name, args=args)
    downloaded = await results.get_file(name)
    if downloaded is None:
        raise FileNotFoundError(f"aria2 downloaded no {name}")
    return downloaded
