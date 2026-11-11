"""canu — Canu de novo assembly of long reads (PacBio CLR, HiFi or Oxford Nanopore).

Exposes :func:`canu`, which assembles one sample's reads in one of Canu's
read modes (``-pacbio``, ``-nanopore``, ``-pacbio-hifi``) for a genome size
such as ``5m``. Further Canu options (``key=value``) go in ``args``.

Canu writes many intermediate files; as upstream, only these are kept, with
the FASTAs gzipped (``gzip -n``):

- ``<prefix>.report``: Canu's report;
- ``<prefix>.contigs.fasta.gz``: the assembly, when Canu assembled contigs;
- ``<prefix>.unassembled.fasta.gz``: reads and tigs left out of contigs;
- ``<prefix>.correctedReads.fasta.gz`` and ``<prefix>.trimmedReads.fasta.gz``:
  when correction or trimming ran (not in HiFi mode);
- ``<prefix>.contigs.layout``, ``<prefix>.contigs.layout.readToTig`` and
  ``<prefix>.contigs.layout.tigInfo``: whichever Canu writes.

``maxThreads`` is the task's CPUs, as upstream. Canu would otherwise size
its jobs to the node's memory, so ``maxMemory`` is the task's memory in GiB;
a ``maxMemory`` in ``args`` comes later and wins.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CANU_IMAGE = "community.wave.seqera.io/library/canu_minimap2_samtools:995ba74386646698"

# Upstream's label asks for 12 CPUs and 72 GB; these defaults fit bacterial
# genomes. Raise them with canu_cmd.override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

MODES = ("-pacbio", "-nanopore", "-pacbio-hifi")

canu_cmd = shell.create(
    name="canu",
    image=CANU_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "mode": str, "genome_size": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        MODE={{inputs.mode}}
        GENOME_SIZE={{inputs.genome_size}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        canu \
            -p "$P" \
            genomeSize="$GENOME_SIZE" \
            maxMemory=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS \
            maxThreads=$TASK_CPUS \
            $MODE "${{READS[@]}}"

        for F in "$P".*.fasta; do gzip -n "$F"; done
        shopt -s nullglob
        for F in "$P".report "$P".contigs.fasta.gz "$P".unassembled.fasta.gz "$P".correctedReads.fasta.gz \
                "$P".trimmedReads.fasta.gz "$P".contigs.layout "$P".contigs.layout.readToTig \
                "$P".contigs.layout.tigInfo; do
            mv "$F" "$OUT"/
        done
    """,
)

env = flyte.TaskEnvironment.from_task("canu", canu_cmd.as_task())


@dataclass
class CanuResult:
    report: File
    assembly: File | None  # <prefix>.contigs.fasta.gz
    contigs: File  # <prefix>.unassembled.fasta.gz, as upstream names it
    corrected_reads: File | None
    corrected_trimmed_reads: File | None
    metadata: File | None  # <prefix>.contigs.layout
    contig_position: File | None  # <prefix>.contigs.layout.readToTig
    contig_info: File | None  # <prefix>.contigs.layout.tigInfo
    results: Dir  # everything above


async def canu(
    reads: File | list[File], mode: str, genome_size: str, prefix: str = "sample", args: str = ""
) -> CanuResult:
    """Assemble ``reads`` in Canu read ``mode`` (e.g. ``-pacbio-hifi``) for ``genome_size`` (e.g. ``5m``)."""
    if mode not in MODES:
        raise ValueError(f"Unrecognised mode to run Canu. Options: {', '.join(MODES)}")
    reads = reads if isinstance(reads, list) else [reads]
    if not reads:
        raise ValueError("canu needs reads")
    results = await canu_cmd(reads=reads, mode=mode, genome_size=genome_size, prefix=prefix, args=args)
    return await read_canu_results(results, prefix)


async def read_canu_results(results: Dir, prefix: str) -> CanuResult:
    """Pick the outputs out of ``canu_cmd``'s results Dir (e.g. from an overridden task)."""

    async def find(suffix: str) -> File | None:
        return await results.get_file(f"{prefix}.{suffix}")

    report, unassembled = await find("report"), await find("unassembled.fasta.gz")
    if report is None or unassembled is None:
        raise FileNotFoundError(f"canu wrote no {prefix}.report or {prefix}.unassembled.fasta.gz")
    return CanuResult(
        report=report,
        assembly=await find("contigs.fasta.gz"),
        contigs=unassembled,
        corrected_reads=await find("correctedReads.fasta.gz"),
        corrected_trimmed_reads=await find("trimmedReads.fasta.gz"),
        metadata=await find("contigs.layout"),
        contig_position=await find("contigs.layout.readToTig"),
        contig_info=await find("contigs.layout.tigInfo"),
        results=results,
    )
