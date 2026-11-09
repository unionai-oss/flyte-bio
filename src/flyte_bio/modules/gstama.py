"""gstama — TAMA (Transcriptome Annotation by Modular Algorithms) for long-read transcripts.

Exposes:

- :func:`gstama_collapse` — ``tama_collapse.py``: collapse one sample's
  mapped transcript reads (SAM/BAM, ``-b BAM`` in ``args`` for BAM) against
  the genome FASTA into a transcript model, as upstream:
  ``<prefix>_collapsed.bed``, ``<prefix>_trans_read.bed``,
  ``<prefix>_local_density_error.txt``, ``<prefix>_polya.txt``,
  ``<prefix>_read.txt``, ``<prefix>_strand_check.txt``,
  ``<prefix>_trans_report.txt`` and, unless variant calling is off,
  ``<prefix>_varcov.txt`` and ``<prefix>_variants.txt``.
- :func:`gstama_merge` — ``tama_merge.py``: merge transcript models (BEDs)
  listed in a TAMA file list (one ``<bed name>\\t<cap flag>\\t<merge
  priorities>\\t<source name>`` row per BED, naming the BEDs by file name),
  with ``-d merge_dup`` as upstream: ``<prefix>.bed``,
  ``<prefix>_gene_report.txt``, ``<prefix>_merge.txt`` and
  ``<prefix>_trans_report.txt``.
- :func:`gstama_polyacleanup` — ``tama_flnc_polya_cleanup.py``: trim
  remaining poly(A) tails from FLNC reads (FASTA), writing gzipped as
  upstream: ``<prefix>.fa.gz``, ``<prefix>_polya_flnc_report.txt.gz`` and
  ``<prefix>_tails.fa.gz``.

Inputs are linked into a working directory under their own names. TAMA is
single-threaded. ``args`` go to the TAMA script.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GSTAMA_IMAGE = "quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0"
# Upstream labels collapse process_medium and the others process_low; TAMA is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

gstama_collapse_cmd = shell.create(
    name="gstama_collapse",
    image=GSTAMA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${BAM[0]}"); ln -s "${BAM[0]}" "$NAME"
        FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"
        tama_collapse.py -s "$NAME" -f "$FA_NAME" -p "$P" $ARGS || exit 1
        shopt -s nullglob
        mv "${P}"_collapsed.bed "${P}"_trans_read.bed "${P}"_*.txt "$OUT"/
    """,
)

gstama_merge_cmd = shell.create(
    name="gstama_merge",
    image=GSTAMA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"beds": list[File], "filelist": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BEDS=({inputs.beds}); LIST=({inputs.filelist})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        for B in "${BEDS[@]}"; do ln -s "$B" "$(basename "$B")"; done
        LIST_NAME=$(basename "${LIST[0]}"); ln -s "${LIST[0]}" "$LIST_NAME"
        tama_merge.py -f "$LIST_NAME" -d merge_dup -p "$P" $ARGS || exit 1
        mv "$P.bed" "${P}_gene_report.txt" "${P}_merge.txt" "${P}_trans_report.txt" "$OUT"/
    """,
)

gstama_polyacleanup_cmd = shell.create(
    name="gstama_polyacleanup",
    image=GSTAMA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$NAME"
        tama_flnc_polya_cleanup.py -f "$NAME" -p "$P" $ARGS || exit 1
        gzip "$P.fa" "${P}_polya_flnc_report.txt" "${P}_tails.fa"
        mv "$P.fa.gz" "${P}_polya_flnc_report.txt.gz" "${P}_tails.fa.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "gstama", gstama_collapse_cmd.as_task(), gstama_merge_cmd.as_task(), gstama_polyacleanup_cmd.as_task()
)


async def required(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class GstamaCollapseResult:
    bed: File  # <prefix>_collapsed.bed: the collapsed transcript models
    bed_trans_reads: File  # <prefix>_trans_read.bed
    local_density_error: File  # <prefix>_local_density_error.txt
    polya: File  # <prefix>_polya.txt
    read: File  # <prefix>_read.txt
    strand_check: File  # <prefix>_strand_check.txt
    trans_report: File  # <prefix>_trans_report.txt
    varcov: File | None  # <prefix>_varcov.txt, with variant calling
    variants: File | None  # <prefix>_variants.txt, with variant calling
    results: Dir  # all of them


async def gstama_collapse(bam: File, fasta: File, prefix: str = "sample", args: str = "") -> GstamaCollapseResult:
    """Collapse mapped transcript reads ``bam`` against genome ``fasta`` (``-b BAM`` in ``args`` for a BAM)."""
    results = await gstama_collapse_cmd(bam=bam, fasta=fasta, prefix=prefix, args=args)
    files = {}
    for key, suffix in (
        ("bed", "collapsed.bed"),
        ("bed_trans_reads", "trans_read.bed"),
        ("local_density_error", "local_density_error.txt"),
        ("polya", "polya.txt"),
        ("read", "read.txt"),
        ("strand_check", "strand_check.txt"),
        ("trans_report", "trans_report.txt"),
    ):
        files[key] = await required(results, f"{prefix}_{suffix}", "tama_collapse.py")
    varcov = await results.get_file(f"{prefix}_varcov.txt")
    variants = await results.get_file(f"{prefix}_variants.txt")
    return GstamaCollapseResult(varcov=varcov, variants=variants, results=results, **files)


@dataclass
class GstamaMergeResult:
    bed: File  # <prefix>.bed: the merged transcript models
    gene_report: File  # <prefix>_gene_report.txt
    merge: File  # <prefix>_merge.txt
    trans_report: File  # <prefix>_trans_report.txt
    results: Dir  # all of them


async def gstama_merge(beds: list[File], filelist: File, prefix: str = "sample", args: str = "") -> GstamaMergeResult:
    """Merge transcript models ``beds`` as ``filelist`` (a TAMA file list naming them) describes."""
    if not beds:
        raise ValueError("tama_merge.py needs at least one BED")
    if f"{prefix}.bed" in {b.name or os.path.basename(b.path) for b in beds}:
        raise ValueError(f"an input and the output are both {prefix}.bed: set another prefix")
    results = await gstama_merge_cmd(beds=beds, filelist=filelist, prefix=prefix, args=args)
    return GstamaMergeResult(
        bed=await required(results, f"{prefix}.bed", "tama_merge.py"),
        gene_report=await required(results, f"{prefix}_gene_report.txt", "tama_merge.py"),
        merge=await required(results, f"{prefix}_merge.txt", "tama_merge.py"),
        trans_report=await required(results, f"{prefix}_trans_report.txt", "tama_merge.py"),
        results=results,
    )


@dataclass
class GstamaPolyaCleanupResult:
    fasta: File  # <prefix>.fa.gz: the reads with poly(A) tails trimmed
    report: File  # <prefix>_polya_flnc_report.txt.gz
    tails: File  # <prefix>_tails.fa.gz: the trimmed tails
    results: Dir  # all of them


async def gstama_polyacleanup(fasta: File, prefix: str = "sample", args: str = "") -> GstamaPolyaCleanupResult:
    """``fasta`` (FLNC reads) with leftover poly(A) tails trimmed: ``<prefix>.fa.gz``."""
    if (fasta.name or os.path.basename(fasta.path)) in {f"{prefix}.fasta", f"{prefix}.fa"}:
        raise ValueError("Input and output names are the same, set another prefix")
    results = await gstama_polyacleanup_cmd(fasta=fasta, prefix=prefix, args=args)
    return GstamaPolyaCleanupResult(
        fasta=await required(results, f"{prefix}.fa.gz", "tama_flnc_polya_cleanup.py"),
        report=await required(results, f"{prefix}_polya_flnc_report.txt.gz", "tama_flnc_polya_cleanup.py"),
        tails=await required(results, f"{prefix}_tails.fa.gz", "tama_flnc_polya_cleanup.py"),
        results=results,
    )
