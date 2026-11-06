"""unzip — extract an archive with 7-Zip.

Wraps ``7za x`` (p7zip), which reads zip and many other formats: as in
upstream's test, a ``.tar.gz`` extracts to its ``.tar``. The extracted files
are written into the output Dir (upstream names that directory after the
sample or the archive's base name).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Upstream's p7zip container.
UNZIP_IMAGE = "quay.io/biocontainers/p7zip:16.02"

# Single-threaded (upstream's process_single); fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


unzip_cmd = shell.create(
    name="unzip",
    image=UNZIP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        A=({inputs.archive})
        ARGS={inputs.args}
        7za x -o{outputs.results}/ $ARGS "${A[0]}"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "unzip",
    unzip_cmd.as_task(),
)


async def unzip(archive: File, args: str = "") -> Dir:
    """``archive`` extracted into a Dir; ``args`` go to ``7za x``."""
    return await unzip_cmd(archive=archive, args=args)
