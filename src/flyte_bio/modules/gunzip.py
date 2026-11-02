"""gunzip — decompress a single gzipped file.

Wraps ``gzip -cd`` (not the ``gunzip`` binary itself: ``gunzip`` writes
the output with the source file's group ownership, ``gzip -cd >``
doesn't). As upstream, the output keeps the input's name without ``.gz``
(``genome.fa.gz`` -> ``genome.fa``): tools that read names or extensions
see what they expect.
"""



import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Seqera wave coreutils image.
GUNZIP_IMAGE = "community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8"

# Single-threaded; fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


gunzip_cmd = shell.create(
    name="gunzip",
    image=GUNZIP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File},
    outputs={"results": Dir},
    script=r"""
        A=({inputs.archive})
        NAME=$(basename "${A[0]}")
        gzip -cd "${A[0]}" > {outputs.results}/"${NAME%.gz}"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "gunzip",
    gunzip_cmd.as_task(),
)


async def gunzip(archive: File) -> File:
    """``archive`` decompressed, named like it without ``.gz``."""
    results = await gunzip_cmd(archive=archive)
    name = archive.path.rstrip("/").rsplit("/", 1)[-1].removesuffix(".gz")
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"gunzip wrote no {name}")
    return out
