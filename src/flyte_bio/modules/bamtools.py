"""bamtools — BAM statistics and splitting with BamTools.

Exposes:

- :func:`bamtools_stats` — ``bamtools stats`` of one BAM: read counts,
  mapping and pairing summaries, to ``<prefix>.bam.stats``.
- :func:`bamtools_split` — merge one or more BAMs (``bamtools merge``) and
  split the result (``bamtools split``) by what ``args`` asks for:
  ``-reference`` (per reference sequence, as ``<prefix>.REF_<name>.bam`` and
  ``<prefix>.REF_unmapped.bam``), ``-mapped``, ``-paired``, ``-tag <TAG>``, ...
  A split can write one BAM per contig, too many files for a Dir, so the BAMs
  come back as one archive, ``<prefix>.split.tar.gz``, with the BAMs at its
  root.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BAMTOOLS_IMAGE = "community.wave.seqera.io/library/bamtools:2.5.2--b804fdcee08dfdd6"
# Upstream labels stats process_single and split process_low; BamTools is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# BamTools may record the input path, so the BAM is read under its own name.
bamtools_stats_cmd = shell.create(
    name="bamtools_stats",
    image=BAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam})
        P={inputs.prefix}
        ARGS={inputs.args}
        bamtools stats -in "${BAM[0]}" $ARGS > {outputs.results}/"$P.bam.stats"
    """,
)

bamtools_split_cmd = shell.create(
    name="bamtools_split",
    image=BAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bams": list[File], "prefix": str, "args": str},
    outputs={"results": Dir},
    script=r"""
        BAMS=({inputs.bams})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        IN=()
        for B in "${BAMS[@]}"; do IN+=(-in "$B"); done
        W=$(mktemp -d); cd "$W"
        mkdir split
        bamtools merge "${IN[@]}" | bamtools split -stub "split/$P" $ARGS
        ls split/*.bam > /dev/null || { echo "bamtools split wrote no BAMs" >&2; exit 1; }
        tar -czf "$OUT/$P.split.tar.gz" -C split .
    """,
)

env = flyte.TaskEnvironment.from_task("bamtools", bamtools_stats_cmd.as_task(), bamtools_split_cmd.as_task())


async def bamtools_stats(bam: File, prefix: str = "sample", args: str = "") -> File:
    """``bamtools stats`` of ``bam``, as ``<prefix>.bam.stats``."""
    results = await bamtools_stats_cmd(bam=bam, prefix=prefix, args=args)
    stats = await results.get_file(f"{prefix}.bam.stats")
    if stats is None:
        raise FileNotFoundError(f"bamtools stats wrote no {prefix}.bam.stats")
    return stats


async def bamtools_split(bams: list[File], args: str, prefix: str = "sample") -> File:
    """Merge ``bams`` and split them as ``args`` asks; returns the archive ``<prefix>.split.tar.gz`` of the BAMs.

    ``args`` must name a split mode (``-reference``, ``-mapped``, ...);
    BamTools has no default one.
    """
    if not bams:
        raise ValueError("bamtools split needs at least one BAM")
    results = await bamtools_split_cmd(bams=bams, prefix=prefix, args=args)
    archive = await results.get_file(f"{prefix}.split.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"bamtools split wrote no {prefix}.split.tar.gz")
    return archive
