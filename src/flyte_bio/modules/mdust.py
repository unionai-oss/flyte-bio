"""mdust — mask low-complexity (dust) regions of DNA sequences.

Exposes :func:`mdust`, which writes ``<prefix>.fasta``: the input FASTA with
low-complexity stretches masked (as ``N``, or as ``args`` says, e.g. ``-m L``
for lowercase; ``-c`` lists their coordinates instead). As upstream, the
prefix defaults to ``<sample>_mdust`` and may not overwrite the input.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MDUST_IMAGE = "quay.io/biocontainers/mdust:2006.10.17--h470a237_1"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

mdust_cmd = shell.create(
    name="mdust",
    image=MDUST_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        mdust "${FA[0]}" $ARGS > {outputs.results}/"$P.fasta" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("mdust", mdust_cmd.as_task())


async def mdust(fasta: File, prefix: str = "sample_mdust", args: str = "") -> File:
    """``fasta`` with its low-complexity regions masked: ``<prefix>.fasta``."""
    if fasta.path.rstrip("/").rsplit("/", 1)[-1] == f"{prefix}.fasta":
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    results = await mdust_cmd(fasta=fasta, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.fasta")
    if found is None:
        raise FileNotFoundError(f"mdust wrote no {prefix}.fasta")
    return found
