"""publish — lay a pipeline's outputs out in one directory, as upstream's ``--outdir``.

Flyte keeps every task output at its own run- and task-specific path. A
pipeline that wants upstream's results tree describes it as a *layout*
(relative path -> File or Dir) and calls :func:`publish`, which copies each
output into place and returns the tree as a :class:`~flyte.io.Dir`, so a
downstream task can take it as an input.

The tree goes to ``outdir`` (e.g. ``s3://my-bucket/rnaseq``; the cluster's
role needs write access) or, without one, to a fresh location in Flyte's
own storage. Copies within one object store are server-side (objects over
S3's 5 GiB copy limit, and copies across stores, are streamed through the
task); nothing is written to local disk.
"""

import asyncio
from collections.abc import Mapping

import flyte
import flyte.storage as storage
from flyte.io import Dir, File

# S3's CopyObject handles objects up to 5 GiB; bigger ones are streamed.
SERVER_SIDE_COPY_MAX_BYTES = 5 * 1024**3
# Copies in flight at once.
MAX_CONCURRENT_COPIES = 16


def protocol(path: str) -> str:
    return path.split("://", 1)[0] if "://" in path else "file"


async def copy_object(src: str, dst: str) -> None:
    """Copy one object, server-side when source and destination share a store."""
    fs = storage.get_underlying_filesystem(path=src)
    if protocol(src) == protocol(dst) and fs.async_impl:
        info = await fs._info(src)
        if (info.get("size") or 0) <= SERVER_SIDE_COPY_MAX_BYTES:
            await fs._cp_file(src, dst)
            return
    await storage.put_stream(storage.get_stream(src), to_path=dst)


def join(base: str, rel: str) -> str:
    return base.rstrip("/") + "/" + rel.strip("/")


def name(f: File) -> str:
    """A File's base name."""
    return f.path.rstrip("/").rsplit("/", 1)[-1]


async def expand(prefix: str, d: Dir) -> dict[str, File]:
    """Layout entries for every file in ``d``, under ``prefix``, keeping their layout inside ``d``.

    For folders that several outputs share (e.g. one Dir per sample whose
    files all belong in ``stringtie/``), where mapping each Dir to the same
    path would clash.
    """
    base = d.path.rstrip("/")
    return {join(prefix, f.path[len(base) :]): f async for f in d.walk()}


async def publish(layout: Mapping[str, File | Dir | None], outdir: str | None = None) -> Dir:
    """Copy each output to ``<outdir>/<relative path>`` and return the tree as a Dir.

    A File lands at its path; a Dir's files land under its path, keeping
    their layout inside the Dir. ``None`` entries (outputs a run didn't
    produce) are skipped. Without ``outdir`` the tree goes to a new location
    in Flyte's storage.
    """
    if outdir is None:
        outdir = flyte.ctx().raw_data_path.get_random_remote_path("published")

    copies: list[tuple[str, str]] = []
    for rel, output in layout.items():
        if output is None:
            continue
        if isinstance(output, Dir):
            copies.extend((f.path, join(outdir, r)) for r, f in (await expand(rel, output)).items())
        else:
            copies.append((output.path, join(outdir, rel)))

    destinations = [dst for _, dst in copies]
    if len(destinations) != len(set(destinations)):
        clashes = sorted({d for d in destinations if destinations.count(d) > 1})
        raise ValueError(f"publish layout puts several outputs at the same path: {clashes}")

    slots = asyncio.Semaphore(MAX_CONCURRENT_COPIES)

    async def copy(src: str, dst: str) -> None:
        async with slots:
            await copy_object(src, dst)

    await asyncio.gather(*(copy(src, dst) for src, dst in copies))
    return Dir.from_existing_remote(outdir)
