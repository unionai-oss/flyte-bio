"""genmap — fast computation of genome mappability.

Exposes :func:`genmap_index`, which builds GenMap's FM index of a FASTA, the
input to ``genmap map``. The index is a Dir of GenMap's ``index.*`` files (the
directory upstream names ``<prefix>``). ``args`` is upstream's ``ext.args``
(e.g. ``-A skew``) and is split on whitespace.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GENMAP_IMAGE = "quay.io/biocontainers/genmap:1.3.0--h1b792b2_1"
# Upstream's process_high; the index of a small genome needs far less. Building
# one for a human genome wants ~ 32 GB.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

# genmap refuses an index directory that already exists, so it builds into a
# fresh one and the files are moved into the output.
genmap_index_cmd = shell.create(
    name="genmap_index",
    image=GENMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d)
        genmap index --fasta-file "${FA[0]}" --index "$W/$P" $ARGS || exit 1
        mv "$W/$P"/* {outputs.index}/
    """,
)

env = flyte.TaskEnvironment.from_task("genmap", genmap_index_cmd.as_task())


async def genmap_index(fasta: File, prefix: str = "sample", args: str = "") -> Dir:
    """GenMap's index of ``fasta``: a Dir of ``index.*`` files."""
    return await genmap_index_cmd(fasta=fasta, prefix=prefix, args=args)
