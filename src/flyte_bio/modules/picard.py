"""picard — Picard tools.

Exposes:

- :func:`picard_markduplicates`, which flags duplicate reads in a
  coordinate-sorted BAM (keeping them, unless told otherwise via ``args``) and
  writes Picard's duplication metrics.
- :func:`picard_collectmultiplemetrics`, which runs Picard's standard set of
  alignment QC collectors over a BAM (or CRAM, with its reference) and
  returns the ``<prefix>.CollectMultipleMetrics.*_metrics`` tables and the
  ``*.pdf`` plots. Upstream pins a newer Picard for it (3.5.0), so it has its
  own environment and ``env`` aggregates both.

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

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

PICARD_IMAGE = "community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6"
COLLECTMULTIPLEMETRICS_IMAGE = "community.wave.seqera.io/library/picard:3.5.0--be20f299a1924f0d"

PICARD_MEMORY_MB = 8192
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory=f"{PICARD_MEMORY_MB}Mi")
# Upstream's process_single.
COLLECTMULTIPLEMETRICS_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

picard_markduplicates_cmd = shell.create(
    name="picard_markduplicates",
    image=PICARD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "fasta": File | None, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
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
            -Xmx{JAVA_HEAP_MB}M \
            MarkDuplicates \
            $ARGS \
            --INPUT "${{BAM[0]}}" \
            --OUTPUT {{outputs.results}}/"$P.bam" \
            "${{REF[@]}}" \
            --METRICS_FILE {{outputs.results}}/"$P.metrics.txt"
    """,
)

# As upstream, the alignment, its index and the reference (with its .fai
# named after it) sit side by side in the working directory.
picard_collectmultiplemetrics_cmd = shell.create(
    name="picard_collectmultiplemetrics",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=COLLECTMULTIPLEMETRICS_RESOURCES,
    inputs={
        "bam": File,
        "bai": File | None,
        "fasta": File | None,
        "fai": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        shopt -s nullglob; BAI=({{inputs.bai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        BAM_LINK=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_LINK"
        if [ ${{#BAI[@]}} -gt 0 ]; then ln -s "${{BAI[0]}}" "$(basename "${{BAI[0]}}")"; fi
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_LINK=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_LINK"
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_LINK.fai"; fi
            REF=(--REFERENCE_SEQUENCE "$FA_LINK")
        fi
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            CollectMultipleMetrics \
            $ARGS \
            --INPUT "$BAM_LINK" \
            --OUTPUT {{outputs.results}}/"$P.CollectMultipleMetrics" \
            "${{REF[@]}}"
    """,
)


# One image per environment: CollectMultipleMetrics' Picard is newer, so `env` aggregates two.
picard_markduplicates_env = flyte.TaskEnvironment.from_task(
    "picard_markduplicates",
    picard_markduplicates_cmd.as_task(),
)
picard_collectmultiplemetrics_env = flyte.TaskEnvironment.from_task(
    "picard_collectmultiplemetrics",
    picard_collectmultiplemetrics_cmd.as_task(),
)
env = flyte.TaskEnvironment(name="picard", depends_on=[picard_markduplicates_env, picard_collectmultiplemetrics_env])


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
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"picard MarkDuplicates wrote no {name}")
        return found

    return MarkDuplicatesResult(bam=await pick(f"{prefix}.bam"), metrics=await pick(f"{prefix}.metrics.txt"))


@dataclass
class CollectMultipleMetricsResult:
    metrics: list[File]  # <prefix>.CollectMultipleMetrics.*_metrics, sorted by name
    pdfs: list[File]  # <prefix>.CollectMultipleMetrics.*.pdf, sorted by name


async def picard_collectmultiplemetrics(
    bam: File,
    prefix: str,
    bai: File | None = None,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> CollectMultipleMetricsResult:
    """Picard's alignment QC metrics for ``bam`` (a CRAM needs ``fasta``)."""
    results = await picard_collectmultiplemetrics_cmd(
        bam=bam,
        bai=bai,
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
    )
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)
    metrics = [f for f in files if f.name.endswith("_metrics")]
    if not metrics:
        raise FileNotFoundError("picard CollectMultipleMetrics wrote no *_metrics files")
    return CollectMultipleMetricsResult(metrics=metrics, pdfs=[f for f in files if f.name.endswith(".pdf")])
