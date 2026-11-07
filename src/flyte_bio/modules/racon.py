"""racon — consensus polishing of an assembly from read-to-assembly overlaps.

Exposes :func:`racon`, which polishes an assembly with the reads and their
PAF overlaps against it, returning ``<prefix>_assembly_consensus.fasta.gz``
(gzipped without a timestamp, as upstream). racon runs on the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RACON_IMAGE = "community.wave.seqera.io/library/racon:1.4.20--5a7230c7590a8e2e"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

racon_cmd = shell.create(
    name="racon",
    image=RACON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "assembly": File, "paf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        READS=({{inputs.reads}})
        ASSEMBLY=({{inputs.assembly}})
        PAF=({{inputs.paf}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        racon -t $TASK_CPUS "${{READS[0]}}" "${{PAF[0]}}" $ARGS "${{ASSEMBLY[0]}}" > "${{P}}_assembly_consensus.fasta" \
            || exit 1
        gzip -n "${{P}}_assembly_consensus.fasta"
        mv "${{P}}_assembly_consensus.fasta.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("racon", racon_cmd.as_task())


async def racon(reads: File, assembly: File, paf: File, prefix: str = "sample", args: str = "") -> File:
    """Polish ``assembly`` with ``reads`` and their overlaps ``paf``: ``<prefix>_assembly_consensus.fasta.gz``."""
    results = await racon_cmd(reads=reads, assembly=assembly, paf=paf, prefix=prefix, args=args)
    name = f"{prefix}_assembly_consensus.fasta.gz"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"racon wrote no {name}")
    return found
