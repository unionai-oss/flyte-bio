"""zip — bundle files into a zip archive with 7-Zip.

Exposes :func:`zip_files`, which archives files into ``<prefix>.zip`` as
upstream does: ``7z a -l`` over the files staged under ``inputs/`` (``-l``
stores what symlinks point to). The archive holds them by name, without the
``inputs/`` directory.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

P7ZIP_IMAGE = "community.wave.seqera.io/library/p7zip:16.02--94efedc6cfea7db9"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

zip_cmd = shell.create(
    name="zip",
    image=P7ZIP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"files": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FILES=({inputs.files})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        mkdir inputs
        for f in "${FILES[@]}"; do ln -s "$f" inputs/; done
        7z a -l $ARGS "$P.zip" ./inputs/* || exit 1
        mv "$P.zip" {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("zip", zip_cmd.as_task())


async def zip_files(files: list[File], prefix: str = "zipped_files", args: str = "") -> File:
    """``<prefix>.zip`` of ``files``."""
    results = await zip_cmd(files=files, prefix=prefix, args=args)
    archive = await results.get_file(f"{prefix}.zip")
    if archive is None:
        raise FileNotFoundError(f"7z wrote no {prefix}.zip")
    return archive
