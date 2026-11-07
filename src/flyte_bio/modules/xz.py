"""xz — decompress a single xz-compressed file.

Exposes :func:`xz_decompress`. As upstream, the output keeps the input's name
without ``.xz`` (``taxa_sqlite.xz`` -> ``taxa_sqlite``), and ``xz`` writes to
stdout (``--stdout``): decompressing in place would trip over the read-only,
symlinked input. ``-T`` is the task's CPUs; xz only decompresses in parallel
when the file has several blocks, as files compressed multi-threaded do.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

XZ_IMAGE = "community.wave.seqera.io/library/xz:5.2.6--5bcc1ddb5aafa435"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi", disk="20Gi")

xz_decompress_cmd = shell.create(
    name="xz_decompress",
    image=XZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        A=({{inputs.archive}})
        ARGS={{inputs.args}}
        NAME=$(basename "${{A[0]}}")
        xz -T $TASK_CPUS --decompress --stdout $ARGS "${{A[0]}}" > {{outputs.results}}/"${{NAME%.xz}}"
    """,
)

env = flyte.TaskEnvironment.from_task("xz", xz_decompress_cmd.as_task())


async def xz_decompress(archive: File, args: str = "") -> File:
    """``archive`` decompressed, named like it without ``.xz``."""
    results = await xz_decompress_cmd(archive=archive, args=args)
    name = archive.path.rstrip("/").rsplit("/", 1)[-1].removesuffix(".xz")
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"xz wrote no {name}")
    return out
