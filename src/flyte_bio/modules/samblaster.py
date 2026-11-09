"""samblaster — mark duplicates and extract discordant and split reads.

Exposes :func:`samblaster`, which runs samblaster over a read-id grouped BAM
(e.g. straight from the aligner) as upstream does:
``samtools view -h <args2> | samblaster <args> | samtools view <args3> -Sb -``,
writing ``<prefix>.bam``. Pass samblaster's options (``-M``,
``--addMateTags``, ``--removeDups``, …) in ``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# samblaster 0.1.26 with samtools 1.19.2.
SAMBLASTER_IMAGE = (
    "quay.io/biocontainers/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:"
    "60ebac4ad9c6530c0d7bf6844f52ec6916e1e0b1-0"
)
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

samblaster_cmd = shell.create(
    name="samblaster",
    image=SAMBLASTER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str, "args2": str, "args3": str},
    defaults={"args": "", "args2": "", "args3": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        ARGS3={{inputs.args3}}
        set -o pipefail
        W=$(mktemp -d); cd "$W"
        samtools view -h $ARGS2 "${{BAM[0]}}" \
            | samblaster $ARGS \
            | samtools view $ARGS3 -Sb - > {{outputs.results}}/"$P.bam" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("samblaster", samblaster_cmd.as_task())


async def samblaster(bam: File, prefix: str, args: str = "", args2: str = "", args3: str = "") -> File:
    """``bam`` with duplicates marked by samblaster: ``<prefix>.bam``.

    ``args`` goes to samblaster, ``args2`` to the ``samtools view`` that reads
    ``bam``, ``args3`` to the one that writes the output.
    """
    if bam.name == f"{prefix}.bam":
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    results = await samblaster_cmd(bam=bam, prefix=prefix, args=args, args2=args2, args3=args3)
    found = await results.get_file(f"{prefix}.bam")
    if found is None:
        raise FileNotFoundError(f"samblaster wrote no {prefix}.bam")
    return found
