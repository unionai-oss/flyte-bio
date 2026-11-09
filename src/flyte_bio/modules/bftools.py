"""bftools — Bio-Formats command-line tools.

Exposes :func:`bftools_showinf`, which reads an image's metadata with
``showinf`` (without reading pixels or checking for upgrades) and writes its
OME-XML as ``<prefix>.xml.gz``. Bio-Formats' Java heap is 80% of the task's
memory.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

BFTOOLS_IMAGE = "community.wave.seqera.io/library/bftools:8.0.0--3d36406c00785003"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# The bftools launchers read their JVM options from BF_FLAGS and the heap size from BF_MAX_MEM.
bftools_showinf_cmd = shell.create(
    name="bftools_showinf",
    image=BFTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        # The image's entrypoint activates its environment; sourced in a function so that
        # the hook's closing `exec "$@"` sees no arguments.
        activate_environment() {{ . /shell-hook.sh; }}
        activate_environment
        IMAGE=({{inputs.image}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        export BF_FLAGS='-XX:+PerfDisableSharedMem' BF_MAX_MEM="{JAVA_HEAP_MB}m"
        showinf -nopix -no-upgrade -omexml-only $ARGS "${{IMAGE[0]}}" | gzip > "$OUT/$P.xml.gz" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("bftools", bftools_showinf_cmd.as_task())


async def bftools_showinf(image: File, prefix: str = "sample", args: str = "") -> File:
    """``image``'s OME-XML metadata: ``<prefix>.xml.gz``."""
    results = await bftools_showinf_cmd(image=image, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.xml.gz")
    if found is None:
        raise FileNotFoundError(f"showinf wrote no {prefix}.xml.gz")
    return found
