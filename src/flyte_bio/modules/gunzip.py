"""gunzip — decompress a single gzipped file.

Wraps ``gzip -cd`` (not the ``gunzip`` binary itself: ``gunzip`` writes
the output with the source file's group ownership, ``gzip -cd >``
doesn't).
"""



import flyte
from flyte.extras import shell
from flyte.io import File

# Seqera wave coreutils image.
GUNZIP_IMAGE = "community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8"

# Single-threaded; fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


gunzip = shell.create(
    name="gunzip",
    image=GUNZIP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File},
    outputs={"out": File},
    script=r"""
        gzip -cd {inputs.archive} > {outputs.out}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "gunzip",
    gunzip.as_task(),
)
