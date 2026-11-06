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
- :func:`samtools_fastq` — convert an alignment's reads back to FASTQ.
- :func:`samtools_collatefastq` — group an alignment's reads by name, then
  convert them to FASTQ (for coordinate-sorted inputs).

The FASTQ conversions write upstream's names: ``<prefix>_1`` and
``<prefix>_2`` for paired reads (``<prefix>_1`` alone for single-end), or
``<prefix>_interleaved`` with ``interleave`` (paired-end only, uncompressed);
``<prefix>_singleton`` for reads whose mate is missing (not written when
interleaving); and ``<prefix>_other`` for reads flagged neither or both
READ1 and READ2. :func:`samtools_fastq` uses ``.fastq.gz`` (``.fastq`` when
interleaved), :func:`samtools_collatefastq` ``.fq.gz`` (``.fq``). As
upstream, ``samtools fastq`` gets ``--threads`` of the task's CPUs minus 1,
and ``collate`` and its ``fastq`` the task's CPUs each.

Each upstream module pins its own samtools build, so the FASTQ conversions
have their own environments; :data:`env` aggregates them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Pinned biocontainer URI (matches the version the RNA-seq pipeline pins).
SAMTOOLS_IMAGE = "community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5"
# The images the upstream fastq and collatefastq modules pin.
SAMTOOLS_FASTQ_IMAGE = "community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4"
SAMTOOLS_COLLATEFASTQ_IMAGE = "community.wave.seqera.io/library/htslib_samtools_xz:ca51bc8a00d53523"

# Most of these are single-threaded I/O over modest inputs; sort is the one
# memory/throughput-sensitive step.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
SORT_RESOURCES = flyte.Resources(cpu=4, memory="12Gi")
# Upstream's process_low for the FASTQ conversions.
FASTQ_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")


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


samtools_fastq_cmd = shell.create(
    name="samtools_fastq",
    image=SAMTOOLS_FASTQ_IMAGE,
    resources=FASTQ_RESOURCES,
    inputs={"alignment": File, "prefix": str, "interleave": bool, "single_end": bool, "args": str},
    defaults={"interleave": False, "single_end": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        THREADS=$(( TASK_CPUS - 1 ))
        if [ {{inputs.interleave}} = true ] && [ {{inputs.single_end}} = false ]; then
            samtools fastq $ARGS --threads $THREADS -0 "$OUT/${{P}}_other.fastq.gz" "${{IN[0]}}" \
                > "$OUT/${{P}}_interleaved.fastq"
        elif [ {{inputs.single_end}} = true ]; then
            samtools fastq $ARGS --threads $THREADS -0 "$OUT/${{P}}_other.fastq.gz" "${{IN[0]}}" \
                -1 "$OUT/${{P}}_1.fastq.gz" -s "$OUT/${{P}}_singleton.fastq.gz"
        else
            samtools fastq $ARGS --threads $THREADS -0 "$OUT/${{P}}_other.fastq.gz" "${{IN[0]}}" \
                -1 "$OUT/${{P}}_1.fastq.gz" -2 "$OUT/${{P}}_2.fastq.gz" -s "$OUT/${{P}}_singleton.fastq.gz"
        fi
    """,
)


# `collate -O` streams name-grouped records into `fastq`; its `.` is the
# prefix for temporary files, so run in a writable work dir. A CRAM needs
# the reference, linked with its .fai into the work dir (htslib may index it).
# `args` goes to collate, `args2` to fastq.
samtools_collatefastq_cmd = shell.create(
    name="samtools_collatefastq",
    image=SAMTOOLS_COLLATEFASTQ_IMAGE,
    resources=FASTQ_RESOURCES,
    inputs={
        "alignment": File,
        "fasta": File | None,
        "fai": File | None,
        "prefix": str,
        "interleave": bool,
        "single_end": bool,
        "args": str,
        "args2": str,
    },
    defaults={"interleave": False, "single_end": False, "args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        set -o pipefail
        W=$(mktemp -d); cd "$W"
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            NAME=$(basename "${{FA[0]}}")
            ln -s "${{FA[0]}}" "$NAME"
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$NAME.fai"; fi
            REF=(--reference "$NAME")
        fi
        if [ {{inputs.interleave}} = true ] && [ {{inputs.single_end}} = false ]; then
            OUTPUT=()
        elif [ {{inputs.single_end}} = true ]; then
            OUTPUT=(-1 "$OUT/${{P}}_1.fq.gz" -s "$OUT/${{P}}_singleton.fq.gz")
        else
            OUTPUT=(-1 "$OUT/${{P}}_1.fq.gz" -2 "$OUT/${{P}}_2.fq.gz" -s "$OUT/${{P}}_singleton.fq.gz")
        fi
        samtools collate $ARGS --threads $TASK_CPUS "${{REF[@]}}" -O "${{IN[0]}}" . \
            | samtools fastq $ARGS2 --threads $TASK_CPUS "${{REF[@]}}" -0 "$OUT/${{P}}_other.fq.gz" \
                "${{OUTPUT[@]}}" > fastq.out
        if [ ${{#OUTPUT[@]}} -eq 0 ]; then mv fastq.out "$OUT/${{P}}_interleaved.fq"; fi
    """,
)


samtools_tools_env = flyte.TaskEnvironment.from_task(
    "samtools_tools",
    samtools_faidx.as_task(),
    samtools_sort.as_task(),
    samtools_view.as_task(),
    samtools_merge_cmd.as_task(),
    samtools_index.as_task(),
    samtools_stats.as_task(),
    samtools_flagstat.as_task(),
    samtools_idxstats.as_task(),
)
samtools_fastq_env = flyte.TaskEnvironment.from_task("samtools_fastq", samtools_fastq_cmd.as_task())
samtools_collatefastq_env = flyte.TaskEnvironment.from_task(
    "samtools_collatefastq", samtools_collatefastq_cmd.as_task()
)
env = flyte.TaskEnvironment(
    name="samtools", depends_on=[samtools_tools_env, samtools_fastq_env, samtools_collatefastq_env]
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


@dataclass
class SamtoolsFastqResult:
    fastq: list[File]  # <prefix>_1 (and <prefix>_2 for paired reads); empty when interleaved
    interleaved: File | None  # with interleave (paired-end only)
    singleton: File | None  # reads whose mate is missing; None when interleaved
    other: File  # reads flagged neither or both READ1 and READ2
    results: Dir  # everything samtools wrote


async def pick_fastq(results: Dir, prefix: str, ext: str, interleaved: bool, single_end: bool) -> SamtoolsFastqResult:
    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"samtools fastq wrote no {name}")
        return found

    if interleaved:
        fastq, singleton = [], None
        interleaved_file = await pick(f"{prefix}_interleaved.{ext}")
    else:
        mates = ["1"] if single_end else ["1", "2"]
        fastq = [await pick(f"{prefix}_{mate}.{ext}.gz") for mate in mates]
        singleton = await pick(f"{prefix}_singleton.{ext}.gz")
        interleaved_file = None
    return SamtoolsFastqResult(
        fastq=fastq,
        interleaved=interleaved_file,
        singleton=singleton,
        other=await pick(f"{prefix}_other.{ext}.gz"),
        results=results,
    )


async def samtools_fastq(
    alignment: File,
    prefix: str,
    interleave: bool = False,
    single_end: bool = False,
    args: str = "",
) -> SamtoolsFastqResult:
    """The reads of a name-grouped ``alignment`` (SAM/BAM/CRAM) as FASTQ, under upstream's names."""
    results = await samtools_fastq_cmd(
        alignment=alignment, prefix=prefix, interleave=interleave, single_end=single_end, args=args
    )
    return await pick_fastq(results, prefix, "fastq", interleave and not single_end, single_end)


async def samtools_collatefastq(
    alignment: File,
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    interleave: bool = False,
    single_end: bool = False,
    args: str = "",
    args2: str = "",
) -> SamtoolsFastqResult:
    """The reads of ``alignment`` grouped by name and written as FASTQ; a CRAM needs ``fasta``.

    ``args`` goes to ``samtools collate``, ``args2`` to ``samtools fastq``.
    """
    results = await samtools_collatefastq_cmd(
        alignment=alignment,
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        interleave=interleave,
        single_end=single_end,
        args=args,
        args2=args2,
    )
    return await pick_fastq(results, prefix, "fq", interleave and not single_end, single_end)
