"""gatk4 — GATK tools.

Exposes:

- :func:`gatk4_createsequencedictionary` — the ``.dict`` sequence dictionary
  of a FASTA.
- :func:`gatk4_intervallisttobed` — convert a Picard interval list to BED.
- :func:`gatk4_markduplicates` — mark duplicates across one or more BAMs
  (e.g. a sample's lanes), writing BAM or, as upstream does, CRAM via
  samtools.

As upstream, each tool gets 80% of the task's memory as Java heap. Inputs
whose names GATK records (the dictionary's ``UR`` field) are linked into the
working directory under their own names first, so outputs match upstream's.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GATK4_IMAGE = "community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3"
# MarkDuplicates' image adds samtools, for the CRAM conversion.
GATK4_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35"

DEFAULT_MEMORY_MB = 6144
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory=f"{DEFAULT_MEMORY_MB}Mi")
JAVA_OPTIONS = f"-Xmx{int(DEFAULT_MEMORY_MB * 0.8)}M -XX:-UsePerfData"
MARKDUPLICATES_MEMORY_MB = 16384
MARKDUPLICATES_RESOURCES = flyte.Resources(cpu=2, memory=f"{MARKDUPLICATES_MEMORY_MB}Mi")
MARKDUPLICATES_JAVA_OPTIONS = f"-Xmx{int(MARKDUPLICATES_MEMORY_MB * 0.8)}M -XX:-UsePerfData"

# GATK writes the dictionary next to the reference, named after it.
gatk4_createsequencedictionary_cmd = shell.create(
    name="gatk4_createsequencedictionary",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        FA=({{inputs.fasta}})
        NAME=$(basename "${{FA[0]}}")
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{FA[0]}}" "$NAME"
        gatk --java-options "{JAVA_OPTIONS}" CreateSequenceDictionary \
            --REFERENCE "$NAME" --URI "$NAME" --TMP_DIR . $ARGS
        mv ./*.dict {{outputs.results}}/
    """,
)

gatk4_intervallisttobed_cmd = shell.create(
    name="gatk4_intervallisttobed",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"intervals": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        IN=({{inputs.intervals}})
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        gatk --java-options "{JAVA_OPTIONS}" IntervalListToBed \
            --INPUT "${{IN[0]}}" --OUTPUT {{outputs.results}}/{{inputs.prefix}}.bed --TMP_DIR . $ARGS
    """,
)


# As upstream: every BAM is an --INPUT, the FASTA (with its index beside it)
# is the reference, and a `.cram` prefix makes MarkDuplicates write BAM that
# samtools then converts to CRAM and indexes (.crai).
gatk4_markduplicates_cmd = shell.create(
    name="gatk4_markduplicates",
    image=GATK4_SAMTOOLS_IMAGE,
    resources=MARKDUPLICATES_RESOURCES,
    inputs={"bams": list[File], "fasta": File | None, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        BAMS=({{inputs.bams}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        INPUTS=(); for b in "${{BAMS[@]}}"; do INPUTS+=(--INPUT "$b"); done
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            ln -s "${{FA[0]}}" genome.fa
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" genome.fa.fai; fi
            REF=(--REFERENCE_SEQUENCE genome.fa)
        fi
        if [[ "$P" == *.cram ]]; then PBAM="${{P%.cram}}.bam"; else PBAM="$OUT/$P"; fi
        gatk --java-options "{MARKDUPLICATES_JAVA_OPTIONS}" MarkDuplicates \
            "${{INPUTS[@]}}" --OUTPUT "$PBAM" --METRICS_FILE "$OUT/$P.metrics" \
            --TMP_DIR . "${{REF[@]}}" $ARGS
        if [[ "$P" == *.cram ]]; then
            samtools view -Ch -T genome.fa -o "$OUT/$P" "$PBAM"
            samtools index "$OUT/$P"
        fi
    """,
)


# One image per environment: MarkDuplicates' differs, so `env` aggregates two.
gatk4_env = flyte.TaskEnvironment.from_task(
    "gatk4_tools",
    gatk4_createsequencedictionary_cmd.as_task(),
    gatk4_intervallisttobed_cmd.as_task(),
)
gatk4_samtools_env = flyte.TaskEnvironment.from_task("gatk4_samtools", gatk4_markduplicates_cmd.as_task())
env = flyte.TaskEnvironment(name="gatk4", depends_on=[gatk4_env, gatk4_samtools_env])


async def only_file(results: Dir, suffix: str, tool: str) -> File:
    found = [f async for f in results.walk() if f.path.endswith(suffix)]
    if len(found) != 1:
        raise FileNotFoundError(f"GATK {tool} wrote {len(found)} *{suffix} files, expected one")
    return found[0]


async def gatk4_createsequencedictionary(fasta: File, args: str = "") -> File:
    """The sequence dictionary of ``fasta``: ``<fasta name without extension>.dict``."""
    results = await gatk4_createsequencedictionary_cmd(fasta=fasta, args=args)
    return await only_file(results, ".dict", "CreateSequenceDictionary")


async def gatk4_intervallisttobed(intervals: File, prefix: str, args: str = "") -> File:
    """``intervals`` (a Picard ``.interval_list``) as ``<prefix>.bed``."""
    results = await gatk4_intervallisttobed_cmd(intervals=intervals, prefix=prefix, args=args)
    return await only_file(results, f"{prefix}.bed", "IntervalListToBed")


@dataclass
class MarkDuplicatesResult:
    alignment: File  # <prefix> (.cram or .bam)
    index: File | None  # .crai for CRAM; BAM output is indexed only if args ask (--CREATE_INDEX true)
    metrics: File  # <prefix>.metrics


async def gatk4_markduplicates(
    bams: list[File],
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> MarkDuplicatesResult:
    """Mark duplicates across ``bams`` into ``<prefix>``: CRAM when it ends in ``.cram`` (needs ``fasta``).

    BAMs are staged under their own names, so they must be distinct.
    """
    if prefix.endswith(".cram") and fasta is None:
        raise ValueError("CRAM output needs the reference fasta")
    results = await gatk4_markduplicates_cmd(bams=bams, fasta=fasta, fai=fai, prefix=prefix, args=args)
    alignment = await results.get_file(prefix)
    metrics = await results.get_file(f"{prefix}.metrics")
    if alignment is None or metrics is None:
        raise FileNotFoundError(f"MarkDuplicates wrote no {prefix} or {prefix}.metrics")
    index = await results.get_file(f"{prefix}.crai") if prefix.endswith(".cram") else None
    return MarkDuplicatesResult(alignment=alignment, index=index, metrics=metrics)
