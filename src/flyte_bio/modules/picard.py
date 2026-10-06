"""picard — Picard tools.

Exposes :func:`picard_markduplicates`, which flags duplicate reads in a
coordinate-sorted BAM (keeping them, unless told otherwise via ``args``) and
writes Picard's duplication metrics.

Outputs are written under their natural names (``<prefix>.bam``,
``<prefix>.metrics.txt``) inside an output Dir. The reference FASTA and its
``.fai`` are staged into separate input dirs, so they're symlinked side by
side first — htsjdk expects the index next to the FASTA, as upstream stages
them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PICARD_IMAGE = "community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6"

PICARD_MEMORY_MB = 8192
# As upstream: give the JVM 80% of the task's memory.
PICARD_XMX_MB = int(PICARD_MEMORY_MB * 0.8)
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory=f"{PICARD_MEMORY_MB}Mi")

# fasta / fai are `list[File]` (0 or 1 item) rather than `File | None` until
# flyteorg/flyte#8118 is deployed: copilot stages a set optional File as a
# bare path, not the per-input dir the shell glob expects.
picard_markduplicates_cmd = shell.create(
    name="picard_markduplicates",
    image=PICARD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "fasta": list[File], "fai": list[File], "prefix": str, "args": str},
    defaults={"fasta": [], "fai": [], "args": ""},
    outputs={"results": Dir},
    script=rf"""
        BAM=({{inputs.bam}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            ln -s "${{FA[0]}}" genome.fa
            [ ${{#FAI[@]}} -gt 0 ] && ln -s "${{FAI[0]}}" genome.fa.fai
            REF=(--REFERENCE_SEQUENCE genome.fa)
        fi
        picard \
            -Xmx{PICARD_XMX_MB}M \
            MarkDuplicates \
            $ARGS \
            --INPUT "${{BAM[0]}}" \
            --OUTPUT {{outputs.results}}/"$P.bam" \
            "${{REF[@]}}" \
            --METRICS_FILE {{outputs.results}}/"$P.metrics.txt"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "picard",
    picard_markduplicates_cmd.as_task(),
)


@dataclass
class MarkDuplicatesResult:
    bam: File
    metrics: File


async def picard_markduplicates(
    bam: File,
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> MarkDuplicatesResult:
    """Mark duplicates in ``bam``; returns ``<prefix>.bam`` and ``<prefix>.metrics.txt``."""
    results = await picard_markduplicates_cmd(
        bam=bam,
        fasta=[fasta] if fasta is not None else [],
        fai=[fai] if fai is not None else [],
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"picard MarkDuplicates wrote no {name}")
        return found

    return MarkDuplicatesResult(bam=await pick(f"{prefix}.bam"), metrics=await pick(f"{prefix}.metrics.txt"))
