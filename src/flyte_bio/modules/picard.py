"""picard — Picard tools.

Exposes:

- :func:`picard_markduplicates`, which flags duplicate reads in a
  coordinate-sorted BAM (keeping them, unless told otherwise via ``args``) and
  writes Picard's duplication metrics.
- :func:`picard_collectmultiplemetrics`, which runs Picard's standard set of
  alignment QC collectors over a BAM (or CRAM, with its reference) and
  returns the ``<prefix>.CollectMultipleMetrics.*_metrics`` tables and the
  ``*.pdf`` plots.
- :func:`picard_addorreplacereadgroups`, which sets the read group of every
  read in a BAM or CRAM, writing ``<prefix>.bam`` (with ``<prefix>.bai`` when
  ``args`` has ``--CREATE_INDEX``) or ``<prefix>.cram``.
- :func:`picard_createsequencedictionary`, which writes a FASTA's sequence
  dictionary, ``<prefix>.dict``.
- :func:`picard_collecthsmetrics`, which writes hybrid-selection metrics,
  ``<prefix>.CollectHsMetrics.coverage_metrics``, for a BAM over bait and
  target intervals (interval lists, or BEDs converted with the reference's
  ``.dict``, as upstream).

Upstream pins a newer Picard (3.5.0) for every tool but MarkDuplicates, so
those share a second environment and ``env`` aggregates both.

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


# As upstream, the output takes the input's extension (bam or cram); a CRAM
# is decoded and written against the reference, its .fai beside it.
picard_addorreplacereadgroups_cmd = shell.create(
    name="picard_addorreplacereadgroups",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "fasta": File | None, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        SUFFIX="${{READS[0]##*.}}"
        W=$(mktemp -d); cd "$W"
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_LINK=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_LINK"
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_LINK.fai"; fi
            REF=(--REFERENCE_SEQUENCE "$FA_LINK")
        fi
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            AddOrReplaceReadGroups \
            $ARGS \
            "${{REF[@]}}" \
            --INPUT "${{READS[0]}}" \
            --OUTPUT {{outputs.results}}/"$P.$SUFFIX"
    """,
)

picard_createsequencedictionary_cmd = shell.create(
    name="picard_createsequencedictionary",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            CreateSequenceDictionary \
            $ARGS \
            --REFERENCE "${{FA[0]}}" \
            --OUTPUT {{outputs.results}}/"$P.dict"
    """,
)

# The metrics header records the command line, so, as upstream stages them,
# the inputs sit in a scratch directory under their basenames (the intervals
# in baits/ and targets/) and are passed as relative paths. BED intervals are
# converted to interval lists beside themselves, named as upstream names them.
picard_collecthsmetrics_cmd = shell.create(
    name="picard_collecthsmetrics",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=COLLECTMULTIPLEMETRICS_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "bait_intervals": File,
        "target_intervals": File,
        "fasta": File | None,
        "fai": File | None,
        "sequence_dict": File | None,
        "gzi": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}})
        BAITS=({{inputs.bait_intervals}}); TARGETS=({{inputs.target_intervals}})
        shopt -s nullglob
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); DICT=({{inputs.sequence_dict}}); GZI=({{inputs.gzi}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        mkdir baits targets
        for f in "${{BAM[0]}}" "${{BAI[0]}}" "${{FA[@]}}" "${{FAI[@]}}" "${{DICT[@]}}" "${{GZI[@]}}"; do
            ln -s "$f" "$(basename "$f")"
        done
        BAIT=baits/$(basename "${{BAITS[0]}}"); ln -s "${{BAITS[0]}}" "$BAIT"
        TARGET=targets/$(basename "${{TARGETS[0]}}"); ln -s "${{TARGETS[0]}}" "$TARGET"
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then REF=(--REFERENCE_SEQUENCE "$(basename "${{FA[0]}}")"); fi
        to_interval_list() {{
            case "$1" in
                *.bed|*.bed.gz)
                    local out
                    out=$(echo "$1" | sed -E 's/.(bed|bed.gz)$/.interval_list/')
                    picard -Xmx{JAVA_HEAP_MB}M BedToIntervalList --INPUT "$1" --OUTPUT "$out" \
                        --SEQUENCE_DICTIONARY "$(basename "${{DICT[0]}}")" --TMP_DIR . >&2
                    echo "$out" ;;
                *) echo "$1" ;;
            esac
        }}
        BAIT=$(to_interval_list "$BAIT")
        TARGET=$(to_interval_list "$TARGET")
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            CollectHsMetrics \
            $ARGS \
            "${{REF[@]}}" \
            --BAIT_INTERVALS "$BAIT" \
            --TARGET_INTERVALS "$TARGET" \
            --INPUT "$(basename "${{BAM[0]}}")" \
            --OUTPUT "$P.CollectHsMetrics.coverage_metrics"
        mv "$P.CollectHsMetrics.coverage_metrics" {{outputs.results}}/
    """,
)


# One image per environment: every tool but MarkDuplicates uses the newer Picard, so `env` aggregates two.
picard_markduplicates_env = flyte.TaskEnvironment.from_task(
    "picard_markduplicates",
    picard_markduplicates_cmd.as_task(),
)
picard_collectmultiplemetrics_env = flyte.TaskEnvironment.from_task(
    "picard_collectmultiplemetrics",
    picard_collectmultiplemetrics_cmd.as_task(),
    picard_addorreplacereadgroups_cmd.as_task(),
    picard_createsequencedictionary_cmd.as_task(),
    picard_collecthsmetrics_cmd.as_task(),
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


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"picard {tool} wrote no {name}")
    return found


@dataclass
class AddOrReplaceReadGroupsResult:
    reads: File  # <prefix>.bam or <prefix>.cram
    bai: File | None  # <prefix>.bai, when args has --CREATE_INDEX (BAM only)


async def picard_addorreplacereadgroups(
    reads: File,
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> AddOrReplaceReadGroupsResult:
    """Set the read group of every read in ``reads`` (a BAM, or a CRAM with ``fasta``).

    The read-group fields (``-LB``, ``-PL``, ``-PU``, ``-SM`` …) go in ``args``.
    """
    suffix = reads.name.rsplit(".", 1)[-1]
    if reads.name == f"{prefix}.{suffix}":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await picard_addorreplacereadgroups_cmd(reads=reads, fasta=fasta, fai=fai, prefix=prefix, args=args)
    return AddOrReplaceReadGroupsResult(
        reads=await pick_file(results, f"{prefix}.{suffix}", "AddOrReplaceReadGroups"),
        bai=await results.get_file(f"{prefix}.bai"),
    )


async def picard_createsequencedictionary(fasta: File, prefix: str, args: str = "") -> File:
    """The sequence dictionary of ``fasta``: ``<prefix>.dict``."""
    results = await picard_createsequencedictionary_cmd(fasta=fasta, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.dict", "CreateSequenceDictionary")


async def picard_collecthsmetrics(
    bam: File,
    bai: File,
    bait_intervals: File,
    target_intervals: File,
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    sequence_dict: File | None = None,
    gzi: File | None = None,
    args: str = "",
) -> File:
    """Hybrid-selection metrics of ``bam``: ``<prefix>.CollectHsMetrics.coverage_metrics``.

    Intervals are Picard interval lists, or BEDs (``.bed``/``.bed.gz``), which need ``sequence_dict``. A
    bgzipped ``fasta`` also needs its ``gzi``.
    """
    results = await picard_collecthsmetrics_cmd(
        bam=bam,
        bai=bai,
        bait_intervals=bait_intervals,
        target_intervals=target_intervals,
        fasta=fasta,
        fai=fai,
        sequence_dict=sequence_dict,
        gzi=gzi,
        prefix=prefix,
        args=args,
    )
    return await pick_file(results, f"{prefix}.CollectHsMetrics.coverage_metrics", "CollectHsMetrics")
