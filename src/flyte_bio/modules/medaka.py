"""medaka — consensus polishing of a draft assembly with Nanopore reads.

Exposes :func:`medaka`, which runs ``medaka_consensus`` on one sample's reads
against its draft assembly and returns the polished consensus as
``<prefix>.fa.gz`` (gzipped without a timestamp, as upstream).

The draft is linked into the working directory first, since medaka writes its
indexes next to it. medaka runs on the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MEDAKA_IMAGE = "quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

medaka_cmd = shell.create(
    name="medaka",
    image=MEDAKA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "assembly": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        DRAFT=({{inputs.assembly}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{DRAFT[0]}}" "$(basename "${{DRAFT[0]}}")"
        medaka_consensus -t $TASK_CPUS $ARGS -i "${{READS[0]}}" -d "$(basename "${{DRAFT[0]}}")" -o ./ || exit 1
        mv consensus.fasta "$P.fa"
        gzip -n "$P.fa"
        mv "$P.fa.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("medaka", medaka_cmd.as_task())


async def medaka(reads: File, assembly: File, prefix: str = "sample", args: str = "") -> File:
    """Polish ``assembly`` with ``reads``: ``<prefix>.fa.gz``."""
    results = await medaka_cmd(reads=reads, assembly=assembly, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.fa.gz")
    if found is None:
        raise FileNotFoundError(f"medaka wrote no {prefix}.fa.gz")
    return found
