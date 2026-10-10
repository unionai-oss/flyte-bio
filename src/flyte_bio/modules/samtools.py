"""samtools — read/write and query high-throughput sequencing alignments.

Exposes the subset of the samtools suite the RNA-seq pipeline leans on:

- :data:`samtools_faidx` — index a FASTA (and emit chromosome sizes).
- :data:`samtools_sort` — coordinate-sort an alignment.
- :data:`samtools_view` — filter or convert an alignment.
- :func:`samtools_merge` — merge sorted BAMs / CRAMs into one.
- :data:`samtools_index` — build a BAM index (``.bai``).
- :data:`samtools_stats` — full alignment statistics.
- :data:`samtools_flagstat` — FLAG-field tallies.
- :data:`samtools_idxstats` — per-reference mapped/unmapped counts.

All tasks share one htslib+samtools biocontainer, so they live in a single
:class:`flyte.TaskEnvironment`.
"""



import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI (matches the version the RNA-seq pipeline pins).
SAMTOOLS_IMAGE = "community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5"

# Most of these are single-threaded I/O over modest inputs; sort is the one
# memory/throughput-sensitive step.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
SORT_RESOURCES = flyte.Resources(cpu=4, memory="12Gi")


# `faidx` writes the index next to its input, so symlink the FASTA into the
# (writable) work dir first — the inputs mount may be read-only. The `.fai`
# columns are derived from sequence content + line wrapping, not the file
# name, so the md5 is reproducible. Chromosome sizes are just its first two
# columns; emit them too (the pipeline needs them for coverage tracks).
samtools_faidx = shell.create(
    name="samtools_faidx",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File},
    outputs={"fai": File, "sizes": File},
    script=r"""
        ln -s {inputs.fasta} genome.fa
        samtools faidx genome.fa
        cut -f 1,2 genome.fa.fai > genome.sizes
        mv genome.fa.fai {outputs.fai}
        mv genome.sizes {outputs.sizes}
    """,
)


# `args` carries extra sort options, e.g. `-n` to sort by read name.
samtools_sort = shell.create(
    name="samtools_sort",
    image=SAMTOOLS_IMAGE,
    resources=SORT_RESOURCES,
    inputs={"bam": File, "args": str},
    defaults={"args": ""},
    outputs={"bam_sorted": File},
    script=r"""
        ARGS={inputs.args}
        samtools sort $ARGS -o {outputs.bam_sorted} {inputs.bam}
    """,
)


# `args` carries the view options, e.g. `-F 0x900 -b` for a BAM of primary
# alignments. Without an output-format flag the output is SAM without header.
# Decoding a CRAM's records needs its reference: pass `fasta` and `fai`.
samtools_view = shell.create(
    name="samtools_view",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "fasta": File | None, "fai": File | None, "args": str},
    defaults={"args": ""},
    outputs={"out": File},
    script=r"""
        BAM=({inputs.bam})
        shopt -s nullglob; FA=({inputs.fasta}); FAI=({inputs.fai}); shopt -u nullglob
        ARGS={inputs.args}
        REF=()
        if [ ${#FA[@]} -gt 0 ]; then
            W=$(mktemp -d)
            ln -s "${FA[0]}" "$W/genome.fa"
            if [ ${#FAI[@]} -gt 0 ]; then ln -s "${FAI[0]}" "$W/genome.fa.fai"; fi
            REF=(--reference "$W/genome.fa")
        fi
        samtools view $ARGS "${REF[@]}" -o {outputs.out} "${BAM[0]}"
    """,
)


# As upstream: the output format follows `--output-fmt` in args, else the
# first input's; a CRAM needs the reference (linked with its .fai into a
# writable place, as htslib may index it).
samtools_merge_cmd = shell.create(
    name="samtools_merge",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"alignments": list[File], "fasta": File | None, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.alignments})
        shopt -s nullglob; FA=({inputs.fasta}); FAI=({inputs.fai}); shopt -u nullglob
        ARGS={inputs.args}
        case " $ARGS " in
            *"--output-fmt cram"*) EXT=cram ;;
            *"--output-fmt sam"*) EXT=sam ;;
            *"--output-fmt bam"*) EXT=bam ;;
            *) EXT="${IN[0]##*.}" ;;
        esac
        W=$(mktemp -d); REF=()
        if [ ${#FA[@]} -gt 0 ]; then
            ln -s "${FA[0]}" "$W/genome.fa"
            if [ ${#FAI[@]} -gt 0 ]; then ln -s "${FAI[0]}" "$W/genome.fa.fai"; fi
            REF=(--reference "$W/genome.fa")
        fi
        samtools merge $ARGS "${REF[@]}" {outputs.results}/{inputs.prefix}.$EXT "${IN[@]}"
    """,
)


# Index reads the (coordinate-sorted) BAM and writes a sidecar. We name the
# output explicitly so it lands in the outputs mount rather than next to the
# read-only input.
samtools_index = shell.create(
    name="samtools_index",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File},
    outputs={"bai": File},
    script=r"""
        samtools index {inputs.bam} {outputs.bai}
    """,
)


# A CRAM needs its reference: pass `fasta` and its `fai`, linked side by side
# (samtools would otherwise try to index the read-only staged FASTA).
samtools_stats = shell.create(
    name="samtools_stats",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "fasta": File | None, "fai": File | None},
    outputs={"stats": File},
    script=r"""
        BAM=({inputs.bam})
        shopt -s nullglob; FA=({inputs.fasta}); FAI=({inputs.fai}); shopt -u nullglob
        REF=()
        if [ ${#FA[@]} -gt 0 ]; then
            W=$(mktemp -d)
            ln -s "${FA[0]}" "$W/genome.fa"
            [ ${#FAI[@]} -gt 0 ] && ln -s "${FAI[0]}" "$W/genome.fa.fai"
            REF=(--reference "$W/genome.fa")
        fi
        samtools stats "${REF[@]}" "${BAM[0]}" > {outputs.stats}
    """,
)


samtools_flagstat = shell.create(
    name="samtools_flagstat",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File},
    outputs={"flagstat": File},
    script=r"""
        samtools flagstat {inputs.bam} > {outputs.flagstat}
    """,
)


# idxstats needs the index present alongside the BAM under a matching name,
# so stage both under one basename before querying.
samtools_idxstats = shell.create(
    name="samtools_idxstats",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File},
    outputs={"idxstats": File},
    script=r"""
        ln -s {inputs.bam} aln.bam
        ln -s {inputs.bai} aln.bam.bai
        samtools idxstats aln.bam > {outputs.idxstats}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "samtools",
    samtools_faidx.as_task(),
    samtools_sort.as_task(),
    samtools_view.as_task(),
    samtools_merge_cmd.as_task(),
    samtools_index.as_task(),
    samtools_stats.as_task(),
    samtools_flagstat.as_task(),
    samtools_idxstats.as_task(),
)


async def samtools_merge(
    alignments: list[File],
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> File:
    """Merge sorted alignments into ``<prefix>.<ext>`` (the inputs' format unless ``--output-fmt`` says).

    Inputs are staged under their own names, so they must be distinct.
    """
    results = await samtools_merge_cmd(alignments=alignments, fasta=fasta, fai=fai, prefix=prefix, args=args)
    merged = [f async for f in results.walk() if f.path.rsplit("/", 1)[-1].startswith(f"{prefix}.")]
    if len(merged) != 1:
        raise FileNotFoundError(f"samtools merge wrote {len(merged)} {prefix}.* files, expected one")
    return merged[0]
