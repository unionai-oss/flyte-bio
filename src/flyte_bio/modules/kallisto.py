"""kallisto — pseudoalignment index for transcript quantification.

Exposes :func:`kallisto_index`, which builds a kallisto index from a
transcript FASTA (optionally gzipped). As upstream, the index is a single file
named ``kallisto``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

KALLISTO_IMAGE = "community.wave.seqera.io/library/kallisto:0.52.0--31c771060d82d25c"
# Upstream labels it process_medium; a human transcriptome needs about 16 GiB.
INDEX_RESOURCES = flyte.Resources(cpu=1, memory="16Gi")

# Upstream passes no thread count, so kallisto indexes with its default (1);
# `-t` in `args` speeds it up.
kallisto_index_cmd = shell.create(
    name="kallisto_index",
    image=KALLISTO_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        kallisto index $ARGS -i kallisto "${FA[0]}"
        mv kallisto "$OUT"/kallisto
    """,
)

env = flyte.TaskEnvironment.from_task("kallisto", kallisto_index_cmd.as_task())


async def kallisto_index(fasta: File, args: str = "") -> File:
    """The kallisto index of ``fasta``: one file named ``kallisto``."""
    results = await kallisto_index_cmd(fasta=fasta, args=args)
    index = await results.get_file("kallisto")
    if index is None:
        raise FileNotFoundError("kallisto index wrote no kallisto")
    return index
