"""unzipfiles — extract an archive with 7-Zip and list its files.

Exposes :func:`unzipfiles`, which runs ``7za x`` (p7zip) on one archive: zip
and many other formats, and, as in upstream's test, a ``.tar.gz`` extracts to
its ``.tar``. Like the ``unzip`` module, but as upstream's ``files`` output it
returns every extracted file, at any depth, in path order, alongside the
extracted Dir (upstream names that directory after the sample or the
archive's base name).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

UNZIPFILES_IMAGE = "quay.io/biocontainers/p7zip:16.02"

# Single-threaded (upstream's process_single); fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

unzipfiles_cmd = shell.create(
    name="unzipfiles",
    image=UNZIPFILES_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        A=({inputs.archive})
        ARGS={inputs.args}
        7za x -o{outputs.results}/ $ARGS "${A[0]}" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("unzipfiles", unzipfiles_cmd.as_task())


@dataclass
class UnzipfilesResult:
    files: list[File]  # every extracted file, in path order
    results: Dir  # the extracted tree


async def unzipfiles(archive: File, args: str = "") -> UnzipfilesResult:
    """``archive`` extracted with ``7za x`` (``args`` go to it): its files and the extracted Dir."""
    results = await unzipfiles_cmd(archive=archive, args=args)
    files = [f async for f in results.walk(recursive=True)]
    return UnzipfilesResult(files=sorted(files, key=lambda f: f.path), results=results)
