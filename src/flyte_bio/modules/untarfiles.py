"""untarfiles — extract a tar archive and list its files.

Exposes :func:`untarfiles`, which runs ``tar -xavf`` on one archive into a
directory named ``prefix`` (by default the archive's name less its last
extension and any ``.tar``, as upstream names it when there's no sample id).
Unlike the ``untar`` module, nothing is stripped: upstream's ``files`` output
is every extracted file, at any depth, which :func:`untarfiles` returns in
path order alongside the extracted Dir. As upstream's glob, hidden files (and
files under hidden directories, such as macOS ``._*`` entries) are left out of
``files``, though they're still in the Dir.

Upstream deprecated this module in favour of ``untar`` (modules commit
bc71bfe, which makes it fail); this is its last working version, modules
commit f46e02b, whose container and tests are unchanged at the pin.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

UNTARFILES_IMAGE = "community.wave.seqera.io/library/grep_sed_tar:db2951cd23a1ffde"

# Single-threaded (upstream's process_single); fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

# `args` goes before the archive, `args2` (e.g. members to extract) after it.
untarfiles_cmd = shell.create(
    name="untarfiles",
    image=UNTARFILES_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File, "prefix": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=r"""
        A=({inputs.archive})
        P={inputs.prefix}
        ARGS={inputs.args}
        ARGS2={inputs.args2}
        mkdir {outputs.results}/"$P"
        tar -C {outputs.results}/"$P" -xavf $ARGS "${A[0]}" $ARGS2 || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("untarfiles", untarfiles_cmd.as_task())


def default_prefix(archive: File) -> str:
    """The archive's name less its last extension and then any ``.tar``, as upstream."""
    name = archive.name
    base = name.rsplit(".", 1)[0] if "." in name else name
    return base.removesuffix(".tar")


def is_hidden(relative_path: str) -> bool:
    """Whether any component of ``relative_path`` is hidden (starts with ``.``), as upstream's glob skips them."""
    return any(part.startswith(".") for part in relative_path.split("/") if part)


@dataclass
class UntarfilesResult:
    files: list[File]  # every extracted file not hidden (as upstream's glob), in path order
    results: Dir  # the extracted tree, <prefix>/


async def untarfiles(archive: File, prefix: str = "", args: str = "", args2: str = "") -> UntarfilesResult:
    """``archive`` extracted with ``tar -xavf`` into ``<prefix>/``: its files and that Dir.

    ``args`` go before the archive and ``args2`` after it. ``prefix`` defaults to the archive's name
    less its extensions (``hello.tar.gz`` → ``hello``).
    """
    prefix = prefix or default_prefix(archive)
    if not prefix or "/" in prefix:
        raise ValueError(f"untarfiles needs a plain directory name as prefix, got {prefix!r}")
    results = await untarfiles_cmd(archive=archive, prefix=prefix, args=args, args2=args2)
    root = f"{results.path.rstrip('/')}/{prefix}"
    extracted = Dir.from_existing_remote(root)
    files = sorted(
        [f async for f in extracted.walk(recursive=True) if not is_hidden(f.path.removeprefix(root))],
        key=lambda f: f.path,
    )
    return UntarfilesResult(files=files, results=extracted)
