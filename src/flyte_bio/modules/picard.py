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
- :func:`picard_mergesamfiles`, which merges BAMs into ``<prefix>.bam``,
  passing them in name order, as upstream (a single BAM is copied as is).
- :func:`picard_sortsam`, which sorts a BAM by ``sort_order``
  (``queryname``, ``coordinate``, ``duplicate`` or ``unsorted``) into
  ``<prefix>.bam``.
- :func:`picard_filtersamreads`, which keeps or drops a BAM's reads by
  ``filter`` (``includeAligned``, ``excludeAligned``, ``includeReadList``,
  ``excludeReadList``, …) into ``<prefix>.bam`` (with ``<prefix>.bai`` when
  ``args`` has ``--CREATE_INDEX true`` and the BAM is coordinate-sorted).
- :func:`picard_collectwgsmetrics`, which writes whole-genome coverage
  metrics, ``<prefix>.CollectWgsMetrics.coverage_metrics``, for a BAM against
  its reference, optionally over an interval list.
- :func:`picard_collectinsertsizemetrics`, which writes a paired-end BAM's
  insert-size metrics, ``<prefix>.txt``, and histogram, ``<prefix>.pdf``.
- :func:`picard_collectrnaseqmetrics`, which writes RNA-seq alignment
  metrics, ``<prefix>.rna_metrics``, against a refFlat annotation (and a
  ``<prefix>`` chart PDF when ``args`` sets ``--CHART_OUTPUT``).
- :func:`picard_bedtointervallist`, which converts a BED to a Picard interval
  list, ``<prefix>.intervallist``, with a sequence dictionary.

Upstream pins a newer Picard (3.5.0) for every tool but MarkDuplicates and
CollectWgsMetrics, so those share a second environment. CollectWgsMetrics
uses a Picard + R image of its own, a third. ``env`` aggregates all three.

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
COLLECTWGSMETRICS_IMAGE = "community.wave.seqera.io/library/picard_r-base:5659104ba2ac5868"

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


# As upstream: the BAMs go to MergeSamFiles in name order (C collation, as
# Groovy sorts paths); a single BAM is the output as is (copied, where
# upstream links it).
picard_mergesamfiles_cmd = shell.create(
    name="picard_mergesamfiles",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bams": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAMS=({{inputs.bams}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        mapfile -t SORTED < <(for f in "${{BAMS[@]}}"; do printf '%s\t%s\n' "$(basename "$f")" "$f"; done \
            | LC_ALL=C sort -t$'\t' -k1,1 | cut -f2)
        if [ ${{#SORTED[@]}} -gt 1 ]; then
            INPUTS=()
            for f in "${{SORTED[@]}}"; do INPUTS+=(--INPUT "$f"); done
            picard \
                -Xmx{JAVA_HEAP_MB}M \
                MergeSamFiles \
                $ARGS \
                "${{INPUTS[@]}}" \
                --OUTPUT "$OUT/$P.bam" || exit 1
        else
            cp "${{SORTED[0]}}" "$OUT/$P.bam"
        fi
    """,
)


picard_sortsam_cmd = shell.create(
    name="picard_sortsam",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "sort_order": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        SORT_ORDER={{inputs.sort_order}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            SortSam \
            --INPUT "${{BAM[0]}}" \
            --OUTPUT {{outputs.results}}/"$P.bam" \
            --SORT_ORDER "$SORT_ORDER" \
            $ARGS || exit 1
    """,
)


# As upstream: the read list is passed only to the read-list filters
# (includeReadList, excludeReadList), and the reference is optional.
picard_filtersamreads_cmd = shell.create(
    name="picard_filtersamreads",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "readlist": File | None,
        "fasta": File | None,
        "filter": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        shopt -s nullglob; READLIST=({{inputs.readlist}}); FA=({{inputs.fasta}}); shopt -u nullglob
        FILTER={{inputs.filter}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        OPTIONAL=()
        if [ ${{#FA[@]}} -gt 0 ]; then OPTIONAL+=(--REFERENCE_SEQUENCE "${{FA[0]}}"); fi
        if [[ "$FILTER" == *ReadList ]]; then OPTIONAL+=(--READ_LIST_FILE "${{READLIST[0]}}"); fi
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            FilterSamReads \
            "${{OPTIONAL[@]}}" \
            --INPUT "${{BAM[0]}}" \
            --OUTPUT {{outputs.results}}/"$P.bam" \
            --FILTER "$FILTER" \
            $ARGS || exit 1
    """,
)


# As upstream, the BAM and its index (when given) sit side by side, as do the
# reference and its .fai.
picard_collectwgsmetrics_cmd = shell.create(
    name="picard_collectwgsmetrics",
    image=COLLECTWGSMETRICS_IMAGE,
    resources=COLLECTMULTIPLEMETRICS_RESOURCES,
    inputs={
        "bam": File,
        "bai": File | None,
        "fasta": File,
        "fai": File | None,
        "intervals": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); FA=({{inputs.fasta}})
        shopt -s nullglob
        BAI=({{inputs.bai}}); FAI=({{inputs.fai}}); INTERVALS=({{inputs.intervals}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for f in "${{BAM[0]}}" "${{BAI[@]}}" "${{FA[0]}}" "${{FAI[@]}}" "${{INTERVALS[@]}}"; do
            ln -s "$f" "$(basename "$f")"
        done
        INTERVAL=()
        if [ ${{#INTERVALS[@]}} -gt 0 ]; then INTERVAL=(--INTERVALS "$(basename "${{INTERVALS[0]}}")"); fi
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            CollectWgsMetrics \
            $ARGS \
            --INPUT "$(basename "${{BAM[0]}}")" \
            --OUTPUT "$P.CollectWgsMetrics.coverage_metrics" \
            --REFERENCE_SEQUENCE "$(basename "${{FA[0]}}")" \
            "${{INTERVAL[@]}}" || exit 1
        mv "$P.CollectWgsMetrics.coverage_metrics" {{outputs.results}}/
    """,
)


picard_collectinsertsizemetrics_cmd = shell.create(
    name="picard_collectinsertsizemetrics",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=COLLECTMULTIPLEMETRICS_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{BAM[0]}}" "$(basename "${{BAM[0]}}")"
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            CollectInsertSizeMetrics \
            --INPUT "$(basename "${{BAM[0]}}")" \
            --OUTPUT "$P.txt" \
            --Histogram_FILE "$P.pdf" \
            $ARGS || exit 1
        mv "$P.txt" "$P.pdf" {{outputs.results}}/
    """,
)


# The reference and the rRNA intervals are optional, as upstream.
picard_collectrnaseqmetrics_cmd = shell.create(
    name="picard_collectrnaseqmetrics",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=COLLECTMULTIPLEMETRICS_RESOURCES,
    inputs={
        "bam": File,
        "ref_flat": File,
        "fasta": File | None,
        "rrna_intervals": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); REF_FLAT=({{inputs.ref_flat}})
        shopt -s nullglob; FA=({{inputs.fasta}}); RRNA=({{inputs.rrna_intervals}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for f in "${{BAM[0]}}" "${{REF_FLAT[0]}}" "${{FA[@]}}" "${{RRNA[@]}}"; do
            ln -s "$f" "$(basename "$f")"
        done
        OPTIONAL=()
        if [ ${{#FA[@]}} -gt 0 ]; then OPTIONAL+=(--REFERENCE_SEQUENCE "$(basename "${{FA[0]}}")"); fi
        if [ ${{#RRNA[@]}} -gt 0 ]; then OPTIONAL+=(--RIBOSOMAL_INTERVALS "$(basename "${{RRNA[0]}}")"); fi
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            CollectRnaSeqMetrics \
            $ARGS \
            "${{OPTIONAL[@]}}" \
            --REF_FLAT "$(basename "${{REF_FLAT[0]}}")" \
            --INPUT "$(basename "${{BAM[0]}}")" \
            --OUTPUT "$P.rna_metrics" || exit 1
        shopt -s nullglob
        mv "$P.rna_metrics" ./*.pdf {{outputs.results}}/
    """,
)


picard_bedtointervallist_cmd = shell.create(
    name="picard_bedtointervallist",
    image=COLLECTMULTIPLEMETRICS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bed": File, "sequence_dict": File, "arguments_file": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BED=({{inputs.bed}}); DICT=({{inputs.sequence_dict}})
        shopt -s nullglob; ARGUMENTS_FILE=({{inputs.arguments_file}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for f in "${{BED[0]}}" "${{DICT[0]}}" "${{ARGUMENTS_FILE[@]}}"; do
            ln -s "$f" "$(basename "$f")"
        done
        ARGS_FILE=()
        if [ ${{#ARGUMENTS_FILE[@]}} -gt 0 ]; then
            ARGS_FILE=(--arguments_file "$(basename "${{ARGUMENTS_FILE[0]}}")")
        fi
        picard \
            -Xmx{JAVA_HEAP_MB}M \
            BedToIntervalList \
            --INPUT "$(basename "${{BED[0]}}")" \
            --OUTPUT "$P.intervallist" \
            --SEQUENCE_DICTIONARY "$(basename "${{DICT[0]}}")" \
            --TMP_DIR . \
            "${{ARGS_FILE[@]}}" \
            $ARGS || exit 1
        mv "$P.intervallist" {{outputs.results}}/
    """,
)


# One image per environment: every tool but MarkDuplicates and CollectWgsMetrics uses the newer Picard, so `env`
# aggregates three.
picard_markduplicates_env = flyte.TaskEnvironment.from_task(
    "picard_markduplicates",
    picard_markduplicates_cmd.as_task(),
)
picard_collectwgsmetrics_env = flyte.TaskEnvironment.from_task(
    "picard_collectwgsmetrics",
    picard_collectwgsmetrics_cmd.as_task(),
)
picard_collectmultiplemetrics_env = flyte.TaskEnvironment.from_task(
    "picard_collectmultiplemetrics",
    picard_collectmultiplemetrics_cmd.as_task(),
    picard_addorreplacereadgroups_cmd.as_task(),
    picard_createsequencedictionary_cmd.as_task(),
    picard_collecthsmetrics_cmd.as_task(),
    picard_mergesamfiles_cmd.as_task(),
    picard_sortsam_cmd.as_task(),
    picard_filtersamreads_cmd.as_task(),
    picard_collectinsertsizemetrics_cmd.as_task(),
    picard_collectrnaseqmetrics_cmd.as_task(),
    picard_bedtointervallist_cmd.as_task(),
)
env = flyte.TaskEnvironment(
    name="picard",
    depends_on=[picard_markduplicates_env, picard_collectmultiplemetrics_env, picard_collectwgsmetrics_env],
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


async def picard_mergesamfiles(bams: list[File], prefix: str, args: str = "") -> File:
    """Merge ``bams`` (in name order, as upstream) into ``<prefix>.bam``; a single BAM is returned as a copy.

    The BAMs are staged under their own names, so they must be distinct.
    """
    names = [f.name for f in bams]
    if not names:
        raise ValueError("picard MergeSamFiles needs at least one BAM")
    if len(set(names)) != len(names):
        raise ValueError(f"picard MergeSamFiles needs distinct BAM names, got {names}")
    if f"{prefix}.bam" in names:
        raise ValueError("an input and the output have the same name; choose another prefix")
    results = await picard_mergesamfiles_cmd(bams=bams, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.bam", "MergeSamFiles")


async def picard_sortsam(bam: File, sort_order: str, prefix: str, args: str = "") -> File:
    """``bam`` sorted by ``sort_order`` (``queryname``, ``coordinate``…): ``<prefix>.bam``."""
    if bam.name == f"{prefix}.bam":
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    results = await picard_sortsam_cmd(bam=bam, sort_order=sort_order, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.bam", "SortSam")


@dataclass
class FilterSamReadsResult:
    bam: File  # <prefix>.bam
    bai: File | None  # <prefix>.bai, when args has --CREATE_INDEX true (coordinate-sorted BAMs only)


async def picard_filtersamreads(
    bam: File,
    filter: str,
    prefix: str,
    readlist: File | None = None,
    fasta: File | None = None,
    args: str = "",
) -> FilterSamReadsResult:
    """``bam`` filtered by ``filter``: ``<prefix>.bam``.

    The read-list filters (``includeReadList``, ``excludeReadList``) need
    ``readlist``, a file of read names.
    """
    if bam.name == f"{prefix}.bam":
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    if filter.endswith("ReadList") and readlist is None:
        raise ValueError(f"picard FilterSamReads --FILTER {filter} needs a read list")
    results = await picard_filtersamreads_cmd(
        bam=bam, readlist=readlist, fasta=fasta, filter=filter, prefix=prefix, args=args
    )
    return FilterSamReadsResult(
        bam=await pick_file(results, f"{prefix}.bam", "FilterSamReads"), bai=await results.get_file(f"{prefix}.bai")
    )


async def picard_collectwgsmetrics(
    bam: File,
    fasta: File,
    prefix: str,
    bai: File | None = None,
    fai: File | None = None,
    intervals: File | None = None,
    args: str = "",
) -> File:
    """Whole-genome coverage metrics of ``bam``: ``<prefix>.CollectWgsMetrics.coverage_metrics``.

    ``intervals``, a Picard interval list, restricts the metrics to its regions.
    """
    results = await picard_collectwgsmetrics_cmd(
        bam=bam, bai=bai, fasta=fasta, fai=fai, intervals=intervals, prefix=prefix, args=args
    )
    return await pick_file(results, f"{prefix}.CollectWgsMetrics.coverage_metrics", "CollectWgsMetrics")


@dataclass
class CollectInsertSizeMetricsResult:
    metrics: File  # <prefix>.txt
    histogram: File  # <prefix>.pdf


async def picard_collectinsertsizemetrics(bam: File, prefix: str, args: str = "") -> CollectInsertSizeMetricsResult:
    """Insert-size metrics of a paired-end ``bam``: ``<prefix>.txt`` and the histogram ``<prefix>.pdf``."""
    results = await picard_collectinsertsizemetrics_cmd(bam=bam, prefix=prefix, args=args)
    return CollectInsertSizeMetricsResult(
        metrics=await pick_file(results, f"{prefix}.txt", "CollectInsertSizeMetrics"),
        histogram=await pick_file(results, f"{prefix}.pdf", "CollectInsertSizeMetrics"),
    )


@dataclass
class CollectRnaSeqMetricsResult:
    metrics: File  # <prefix>.rna_metrics
    pdfs: list[File]  # the chart, when args sets --CHART_OUTPUT


async def picard_collectrnaseqmetrics(
    bam: File,
    ref_flat: File,
    prefix: str,
    fasta: File | None = None,
    rrna_intervals: File | None = None,
    args: str = "",
) -> CollectRnaSeqMetricsResult:
    """RNA-seq alignment metrics of ``bam`` against the ``ref_flat`` annotation: ``<prefix>.rna_metrics``.

    The strand (``--STRAND_SPECIFICITY``) goes in ``args``; ``rrna_intervals`` is a Picard interval list of
    the rRNA genes.
    """
    results = await picard_collectrnaseqmetrics_cmd(
        bam=bam, ref_flat=ref_flat, fasta=fasta, rrna_intervals=rrna_intervals, prefix=prefix, args=args
    )
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)
    return CollectRnaSeqMetricsResult(
        metrics=await pick_file(results, f"{prefix}.rna_metrics", "CollectRnaSeqMetrics"),
        pdfs=[f for f in files if f.name.endswith(".pdf")],
    )


async def picard_bedtointervallist(
    bed: File, sequence_dict: File, prefix: str, arguments_file: File | None = None, args: str = ""
) -> File:
    """``bed`` as a Picard interval list with ``sequence_dict``'s header: ``<prefix>.intervallist``.

    ``arguments_file`` is a file of further Picard arguments (``--arguments_file``).
    """
    results = await picard_bedtointervallist_cmd(
        bed=bed, sequence_dict=sequence_dict, arguments_file=arguments_file, prefix=prefix, args=args
    )
    return await pick_file(results, f"{prefix}.intervallist", "BedToIntervalList")
