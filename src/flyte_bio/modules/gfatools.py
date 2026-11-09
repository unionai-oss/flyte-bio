"""gfatools — tools for GFA assembly graphs.

Exposes :func:`gfatools_gfa2fa`, which writes the segment sequences of a GFA
(e.g. hifiasm's contig graphs) as FASTA: ``<prefix>.fasta.gz``, as upstream.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GFATOOLS_IMAGE = "quay.io/biocontainers/gfatools:0.5--h577a1d6_5"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

gfatools_gfa2fa_cmd = shell.create(
    name="gfatools_gfa2fa",
    image=GFATOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gfa": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        GFA=({inputs.gfa})
        ARGS={inputs.args}
        gfatools gfa2fa $ARGS "${GFA[0]}" | gzip -c > {outputs.results}/{inputs.prefix}.fasta.gz || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("gfatools", gfatools_gfa2fa_cmd.as_task())


async def gfatools_gfa2fa(gfa: File, prefix: str = "sample", args: str = "") -> File:
    """The segment sequences of ``gfa``: ``<prefix>.fasta.gz``."""
    results = await gfatools_gfa2fa_cmd(gfa=gfa, prefix=prefix, args=args)
    fasta = await results.get_file(f"{prefix}.fasta.gz")
    if fasta is None:
        raise FileNotFoundError(f"gfatools gfa2fa wrote no {prefix}.fasta.gz")
    return fasta
