"""pigz — parallel gzip compression and decompression.

Exposes:

- :func:`pigz_compress` — gzip a file to ``<name>.gz``.
- :func:`pigz_uncompress` — decompress a gzipped file to its name without
  ``.gz`` (``test_1.fastq.gz`` -> ``test_1.fastq``).

Both run on the task's CPUs (``pigz -p``).
"""

import os

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PIGZ_IMAGE = "quay.io/biocontainers/pigz:2.8"

# Upstream's process_low gives pigz 2 CPUs; it needs little memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="2Gi")

# As upstream, --stdout: pigz otherwise skips a symlinked input.
pigz_compress_cmd = shell.create(
    name="pigz_compress",
    image=PIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"raw_file": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.raw_file}})
        ARGS={{inputs.args}}
        NAME=$(basename "${{IN[0]}}")
        pigz --processes $TASK_CPUS --stdout --force $ARGS "${{IN[0]}}" > {{outputs.results}}/"$NAME".gz
    """,
)

# As upstream, `unpigz -fk` beside the input: here a symlink to it in the
# output directory (inputs are read-only), removed afterwards.
pigz_uncompress_cmd = shell.create(
    name="pigz_uncompress",
    image=PIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.archive}})
        ARGS={{inputs.args}}
        NAME=$(basename "${{IN[0]}}")
        cd {{outputs.results}}
        ln -s "${{IN[0]}}" "$NAME"
        unpigz -p $TASK_CPUS -fk $ARGS "$NAME"
        rm "$NAME"
    """,
)


env = flyte.TaskEnvironment.from_task("pigz", pigz_compress_cmd.as_task(), pigz_uncompress_cmd.as_task())


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path.rstrip("/"))


async def pigz_compress(raw_file: File, args: str = "") -> File:
    """``raw_file`` gzipped, as ``<name>.gz``. ``args`` go to pigz."""
    name = f"{file_name(raw_file)}.gz"
    results = await pigz_compress_cmd(raw_file=raw_file, args=args)
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"pigz wrote no {name}")
    return out


async def pigz_uncompress(archive: File, args: str = "") -> File:
    """``archive`` decompressed, named like it without ``.gz``. ``args`` go to unpigz."""
    name = file_name(archive).removesuffix(".gz")
    results = await pigz_uncompress_cmd(archive=archive, args=args)
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"unpigz wrote no {name}")
    return out
