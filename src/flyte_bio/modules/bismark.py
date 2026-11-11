"""bismark — bisulfite read alignment and methylation calling with Bismark.

Exposes:

- :func:`bismark_genomepreparation` — the bisulfite index of a FASTA
  (``bismark_genome_preparation``): a Dir holding the FASTA and its
  ``Bisulfite_Genome/{CT,GA}_conversion/`` indexes, built for Bowtie 2 by
  default (``--hisat2``, ``--minimap2`` in ``args`` pick another aligner).
- :func:`bismark_align` — align one sample's reads to that index
  (``bismark --bam``): ``<reads>_bismark_<aligner>[_pe].bam``, its report
  ``*_report.txt`` and, with ``--unmapped`` or ``--ambiguous``, ``*.fq.gz``.
- :func:`bismark_deduplicate` — remove duplicate alignments
  (``deduplicate_bismark``): ``<bam name>.deduplicated.bam`` and
  ``<bam name>.deduplication_report.txt``.
- :func:`bismark_methylationextractor` — methylation calls
  (``bismark_methylation_extractor --bedGraph --counts --gzip --report``):
  per-context, per-strand calls ``<context>_<strand>_<bam name>.txt.gz``, the
  ``<bam name>.bedGraph.gz`` and ``<bam name>.bismark.cov.gz`` summaries,
  ``<bam name>_splitting_report.txt`` and ``<bam name>.M-bias.txt``.
- :func:`bismark_coverage2cytosine` — a genome-wide cytosine report from a
  coverage file (``coverage2cytosine --gzip``):
  ``<prefix>.CpG_report.txt.gz`` (``*report.txt.gz``),
  ``<prefix>.cytosine_context_summary.txt`` and, with options such as
  ``--merge_CpG``, ``*.cov.gz``.
- :func:`bismark_report` — one sample's HTML report from its alignment,
  deduplication, splitting and M-bias reports (``bismark2report``):
  ``<reads>_bismark_<aligner>_{SE,PE}_report.html``.
- :func:`bismark_summary` — the summary across samples (``bismark2summary``):
  ``bismark_summary_report.{txt,html}``.

As upstream:

- The index is linked into the working directory as ``BismarkIndex``, the
  directory upstream's genome preparation writes, so options in ``args`` can
  name it (``--genome_folder BismarkIndex`` for the extractor's
  ``--cytosine_report``).
- :func:`bismark_align` runs one aligner instance per strand (2, or 4 with
  ``--non_directional``; 1 with ``--combined_index``, 2 with
  ``--combined_index_parallel``) and gives each ``-p <task CPUs / instances>``
  threads when that is at least 2. With ``--minimap2`` / ``--rammap`` it
  instead splits the reads into ``--multicore`` chunks: one per 3 CPUs and
  13 GiB (5 CPUs and 18 GiB with ``--non_directional``). Options already in
  ``args`` (``--multicore``, ``--parallel``, ``-p``) win.
- The methylation extractor gets ``--multicore <task CPUs / 3>`` with at least
  6 CPUs, and ``--buffer_size <memory - 2>G`` with more than 6 GiB.
- Inputs keep their names, from which Bismark names its outputs.

Upstream's alignment, extractor and cytosine-report steps also stage the
FASTA, which their commands don't read (the index holds a copy); the wrappers
don't take it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BISMARK_IMAGE = "community.wave.seqera.io/library/bismark:3.1.0--9557d6ab108a83e4"

# Upstream labels genome preparation, alignment, deduplication and extraction
# process_high, the rest process_low. These defaults fit small genomes; raise
# them with <cmd>.override(resources=...) for larger ones. The 4 CPUs
# match upstream's tests (alignment records its per-instance -p in its report).
HIGH_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
LOW_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

bismark_genomepreparation_cmd = shell.create(
    name="bismark_genomepreparation",
    image=BISMARK_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        ARGS={{inputs.args}}
        OUT={{outputs.index}}
        cp "${{FA[0]}}" "$OUT/$(basename "${{FA[0]}}")"
        bismark_genome_preparation $ARGS "$OUT"
    """,
)

bismark_align_cmd = shell.create(
    name="bismark_align",
    image=BISMARK_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "index": Dir, "prefix": str, "args": str},
    defaults={"prefix": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        INDEX_DIR={{inputs.index}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s "$INDEX_DIR" BismarkIndex
        if [ -n "$P" ]; then ARGS="$ARGS --prefix $P"; fi
        ln -s "${{R1[0]}}" "$(basename "${{R1[0]}}")"
        if [ ${{#R2[@]}} -gt 0 ]; then
            ln -s "${{R2[0]}}" "$(basename "${{R2[0]}}")"
            READS=(-1 "$(basename "${{R1[0]}}")" -2 "$(basename "${{R2[0]}}")")
        else
            READS=("$(basename "${{R1[0]}}")")
        fi

        has() {{ [[ "$ARGS" == *"$1"* ]]; }}
        if has --minimap2 || has --mm2 || has --rammap || has --ram; then
            if ! has --multicore; then
                CPU_PER_CORE=3; MB_PER_CORE=$(( 13 * 1024 ))
                if has --non_directional; then CPU_PER_CORE=5; MB_PER_CORE=$(( 18 * 1024 )); fi
                CORES=$(( TASK_CPUS / CPU_PER_CORE ))
                MEMORY_CORES=$(( TASK_MEMORY_MB / MB_PER_CORE ))
                if [ $MEMORY_CORES -lt $CORES ]; then CORES=$MEMORY_CORES; fi
                if [ $CORES -gt 1 ]; then ARGS="$ARGS --multicore $CORES"; fi
            fi
        elif ! has --multicore && ! has --parallel && ! has "-p "; then
            INSTANCES=2
            if has --combined_index; then
                INSTANCES=1
                if has --non_directional && has --combined_index_parallel; then INSTANCES=2; fi
            elif has --non_directional; then
                INSTANCES=4
            fi
            THREADS=$(( TASK_CPUS / INSTANCES ))
            if [ $THREADS -ge 2 ]; then ARGS="$ARGS -p $THREADS"; fi
        fi

        bismark "${{READS[@]}}" --genome BismarkIndex --bam $ARGS
        shopt -s nullglob
        mv ./*bam ./*report.txt ./*fq.gz "$OUT"/
    """,
)

bismark_deduplicate_cmd = shell.create(
    name="bismark_deduplicate",
    image=BISMARK_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"bam": File, "single_end": bool, "args": str},
    defaults={"single_end": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        ARGS={{inputs.args}}
        SEQTYPE=-p
        if [ {{inputs.single_end}} = true ]; then SEQTYPE=-s; fi
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{BAM[0]}}")
        ln -s "${{BAM[0]}}" "$NAME"
        deduplicate_bismark $ARGS $SEQTYPE --bam "$NAME"
        mv ./*.deduplicated.bam ./*.deduplication_report.txt {{outputs.results}}/
    """,
)

bismark_methylationextractor_cmd = shell.create(
    name="bismark_methylationextractor",
    image=BISMARK_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"bam": File, "index": Dir | None, "single_end": bool, "args": str},
    defaults={"single_end": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        INDEX_DIR={{inputs.index}}
        ARGS={{inputs.args}}
        SEQTYPE=-p
        if [ {{inputs.single_end}} = true ]; then SEQTYPE=-s; fi
        W=$(mktemp -d); cd "$W"
        if [ -n "$INDEX_DIR" ]; then ln -s "$INDEX_DIR" BismarkIndex; fi
        mkdir input
        NAME=$(basename "${{BAM[0]}}")
        ln -s "${{BAM[0]}}" "input/$NAME"
        if [[ "$ARGS" != *--multicore* ]] && [ $TASK_CPUS -ge 6 ]; then
            ARGS="$ARGS --multicore $(( TASK_CPUS / 3 ))"
        fi
        MEMORY_GB=$(( TASK_MEMORY_MB / 1024 ))
        if [[ "$ARGS" != *--buffer_size* ]] && [ $MEMORY_GB -gt 6 ]; then
            ARGS="$ARGS --buffer_size $(( MEMORY_GB - 2 ))G"
        fi
        bismark_methylation_extractor "input/$NAME" --bedGraph --counts --gzip --report $SEQTYPE $ARGS
        shopt -s nullglob
        mv ./*.bedGraph.gz ./*.txt.gz ./*.cov.gz ./*_splitting_report.txt ./*.M-bias.txt {{outputs.results}}/
    """,
)

bismark_coverage2cytosine_cmd = shell.create(
    name="bismark_coverage2cytosine",
    image=BISMARK_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"coverage": File, "index": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        COV=({{inputs.coverage}})
        INDEX_DIR={{inputs.index}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "$INDEX_DIR" BismarkIndex
        mkdir input
        NAME=$(basename "${{COV[0]}}")
        ln -s "${{COV[0]}}" "input/$NAME"
        coverage2cytosine "input/$NAME" --genome BismarkIndex --output "$P" --gzip $ARGS
        shopt -s nullglob
        mv ./*.cov.gz ./*report.txt.gz ./*cytosine_context_summary.txt {{outputs.results}}/
    """,
)

bismark_report_cmd = shell.create(
    name="bismark_report",
    image=BISMARK_IMAGE,
    resources=LOW_RESOURCES,
    inputs={
        "align_report": File,
        "dedup_report": File | None,
        "splitting_report": File | None,
        "mbias": File | None,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        REPORTS=({{inputs.align_report}} {{inputs.dedup_report}} {{inputs.splitting_report}} {{inputs.mbias}})
        shopt -u nullglob
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for F in "${{REPORTS[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        bismark2report $ARGS
        shopt -s nullglob
        for F in ./*report.html ./*report.txt; do
            if [ ! -L "$F" ]; then mv "$F" {{outputs.results}}/; fi
        done
    """,
)

bismark_summary_cmd = shell.create(
    name="bismark_summary",
    image=BISMARK_IMAGE,
    resources=LOW_RESOURCES,
    inputs={
        "bam_names": str,
        "align_reports": list[File],
        "dedup_reports": list[File],
        "splitting_reports": list[File],
        "mbias": list[File],
    },
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        REPORTS=({{inputs.align_reports}} {{inputs.dedup_reports}} {{inputs.splitting_reports}} {{inputs.mbias}})
        shopt -u nullglob
        BAMS={{inputs.bam_names}}
        W=$(mktemp -d); cd "$W"
        for F in "${{REPORTS[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        bismark2summary $BAMS
        mv bismark_summary_report.txt bismark_summary_report.html {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "bismark",
    bismark_genomepreparation_cmd.as_task(),
    bismark_align_cmd.as_task(),
    bismark_deduplicate_cmd.as_task(),
    bismark_methylationextractor_cmd.as_task(),
    bismark_coverage2cytosine_cmd.as_task(),
    bismark_report_cmd.as_task(),
    bismark_summary_cmd.as_task(),
)


async def files_ending(results: Dir, suffix: str) -> list[File]:
    """The files directly in ``results`` whose names end with ``suffix``, sorted by name."""
    found = [f async for f in results.walk() if f.name.endswith(suffix)]
    return sorted(found, key=lambda f: f.name)


async def one_file(results: Dir, suffix: str, tool: str) -> File:
    found = await files_ending(results, suffix)
    if len(found) != 1:
        raise FileNotFoundError(f"{tool} wrote {len(found)} *{suffix} files, expected one")
    return found[0]


async def bismark_genomepreparation(fasta: File, args: str = "") -> Dir:
    """The bisulfite index of ``fasta``: the FASTA and its ``Bisulfite_Genome/`` (Bowtie 2 unless ``args`` say)."""
    return await bismark_genomepreparation_cmd(fasta=fasta, args=args)


@dataclass
class BismarkAlignResult:
    bam: File
    report: File  # the alignment report, *_report.txt
    unmapped: list[File]  # *.fq.gz, with --unmapped or --ambiguous
    results: Dir  # everything above


async def bismark_align(reads: list[File], index: Dir, prefix: str = "", args: str = "") -> BismarkAlignResult:
    """Align one sample's reads (one file, or an R1/R2 pair) to a :func:`bismark_genomepreparation` index.

    ``prefix`` (``--prefix``) is prepended to the output names; by default they
    follow the reads' names. Pick the aligner in ``args`` (``--hisat2``,
    ``--minimap2``), matching the index.
    """
    if len(reads) not in (1, 2):
        raise ValueError(f"bismark_align takes one or two read files, got {len(reads)}")
    results = await bismark_align_cmd(
        reads_1=reads[0], reads_2=reads[1] if len(reads) == 2 else None, index=index, prefix=prefix, args=args
    )
    return BismarkAlignResult(
        bam=await one_file(results, ".bam", "bismark"),
        report=await one_file(results, "_report.txt", "bismark"),
        unmapped=await files_ending(results, ".fq.gz"),
        results=results,
    )


@dataclass
class BismarkDeduplicateResult:
    bam: File  # <bam name>.deduplicated.bam
    report: File  # <bam name>.deduplication_report.txt
    results: Dir


async def bismark_deduplicate(bam: File, single_end: bool = False, args: str = "") -> BismarkDeduplicateResult:
    """Remove duplicate alignments from a Bismark BAM."""
    results = await bismark_deduplicate_cmd(bam=bam, single_end=single_end, args=args)
    return BismarkDeduplicateResult(
        bam=await one_file(results, ".deduplicated.bam", "deduplicate_bismark"),
        report=await one_file(results, ".deduplication_report.txt", "deduplicate_bismark"),
        results=results,
    )


@dataclass
class BismarkMethylationExtractorResult:
    bedgraph: File  # <bam name>.bedGraph.gz
    methylation_calls: list[File]  # <context>_<strand>_<bam name>.txt.gz (and any other *.txt.gz)
    coverage: File  # <bam name>.bismark.cov.gz
    report: File  # <bam name>_splitting_report.txt
    mbias: File  # <bam name>.M-bias.txt
    results: Dir


async def bismark_methylationextractor(
    bam: File, index: Dir | None = None, single_end: bool = False, args: str = ""
) -> BismarkMethylationExtractorResult:
    """Extract methylation calls from a Bismark BAM.

    ``index`` is linked as ``BismarkIndex``, for options in ``args`` that read
    the genome (``--cytosine_report --genome_folder BismarkIndex``).
    """
    results = await bismark_methylationextractor_cmd(bam=bam, index=index, single_end=single_end, args=args)
    tool = "bismark_methylation_extractor"
    return BismarkMethylationExtractorResult(
        bedgraph=await one_file(results, ".bedGraph.gz", tool),
        methylation_calls=await files_ending(results, ".txt.gz"),
        coverage=await one_file(results, ".cov.gz", tool),
        report=await one_file(results, "_splitting_report.txt", tool),
        mbias=await one_file(results, ".M-bias.txt", tool),
        results=results,
    )


@dataclass
class BismarkCoverage2CytosineResult:
    report: File  # <prefix>.CpG_report.txt.gz (*report.txt.gz)
    summary: File  # <prefix>.cytosine_context_summary.txt
    coverage: list[File]  # *.cov.gz, with options such as --merge_CpG
    results: Dir


async def bismark_coverage2cytosine(
    coverage: File, index: Dir, prefix: str, args: str = ""
) -> BismarkCoverage2CytosineResult:
    """A genome-wide cytosine report from a Bismark coverage file and the index's genome."""
    results = await bismark_coverage2cytosine_cmd(coverage=coverage, index=index, prefix=prefix, args=args)
    return BismarkCoverage2CytosineResult(
        report=await one_file(results, "report.txt.gz", "coverage2cytosine"),
        summary=await one_file(results, "cytosine_context_summary.txt", "coverage2cytosine"),
        coverage=await files_ending(results, ".cov.gz"),
        results=results,
    )


async def bismark_report(
    align_report: File,
    dedup_report: File | None = None,
    splitting_report: File | None = None,
    mbias: File | None = None,
    args: str = "",
) -> Dir:
    """One sample's ``bismark2report``: a Dir of ``*report.html`` (and ``*report.txt``, if it writes one)."""
    return await bismark_report_cmd(
        align_report=align_report, dedup_report=dedup_report, splitting_report=splitting_report, mbias=mbias, args=args
    )


@dataclass
class BismarkSummaryResult:
    html: File  # bismark_summary_report.html
    txt: File  # bismark_summary_report.txt
    results: Dir


async def bismark_summary(
    bams: list[File],
    align_reports: list[File],
    dedup_reports: list[File],
    splitting_reports: list[File],
    mbias: list[File],
) -> BismarkSummaryResult:
    """The summary across samples. Only the BAMs' names are used: the reports are found from them."""
    results = await bismark_summary_cmd(
        bam_names=" ".join(b.name for b in bams),
        align_reports=align_reports,
        dedup_reports=dedup_reports,
        splitting_reports=splitting_reports,
        mbias=mbias,
    )
    return BismarkSummaryResult(
        html=await one_file(results, "bismark_summary_report.html", "bismark2summary"),
        txt=await one_file(results, "bismark_summary_report.txt", "bismark2summary"),
        results=results,
    )
