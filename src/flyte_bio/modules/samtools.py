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
- :func:`samtools_convert` — convert a BAM to CRAM or a CRAM to BAM against
  its reference, and index the result.
- :func:`samtools_dict` — a FASTA's sequence dictionary, ``<fasta name>.dict``.
- :func:`samtools_calmd` — recalculate an alignment's MD/NM tags against its
  reference, as ``<prefix>.bam`` (SAM text unless ``args`` has ``-b``).
- :func:`samtools_fixmate` — fill in mate coordinates and flags of a
  name-grouped alignment, as ``<prefix>.<ext>``.
- :func:`samtools_markdup` — mark duplicates in a coordinate-sorted
  alignment, as ``<prefix>.<ext>``.
- :func:`samtools_import` — FASTQ reads (single-end, a pair, or one
  interleaved file) as an unaligned ``<prefix>.<ext>``.
- :func:`samtools_trimheader` — an indexed BAM whose header keeps only the
  references with reads (mapped or unmapped, from ``samtools idxstats``), as
  ``<prefix>.bam``. ``args`` goes to the final ``samtools view``.
- :func:`samtools_bam2fq` — an alignment's reads as gzipped FASTQ, split
  into ``<prefix>_1``, ``<prefix>_2``, ``<prefix>_other`` and
  ``<prefix>_singleton`` (``.fq.gz``), or interleaved into
  ``<prefix>_interleaved.fq.gz`` (bgzipped). As upstream, ``-@`` is the
  task's CPUs.
- :func:`samtools_mpileup` — the pileup of an alignment (optionally against
  its reference, within intervals or, with ``-r`` in ``args``, a region of an
  indexed alignment), as ``<prefix>.mpileup.gz``.

fixmate, markdup and import write SAM, BAM or CRAM following ``--output-fmt`` in
``args`` (BAM by default). As upstream, calmd and markdup get ``-@`` of the
task's CPUs, fixmate ``--threads`` of its CPUs minus 1.

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
have their own environments (convert, dict, calmd, fixmate, markdup and
import share fastq's samtools 1.25, bam2fq and mpileup share collatefastq's
image, and trimheader has its own image with gawk); :data:`env` aggregates
them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Pinned biocontainer URI (matches the version the RNA-seq pipeline pins).
SAMTOOLS_IMAGE = "community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5"
# The images the upstream fastq and collatefastq modules pin; convert and
# dict pin fastq's.
SAMTOOLS_FASTQ_IMAGE = "community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4"
SAMTOOLS_COLLATEFASTQ_IMAGE = "community.wave.seqera.io/library/htslib_samtools_xz:ca51bc8a00d53523"
# The image the upstream trimheader module pins (samtools 1.25 with gawk).
SAMTOOLS_TRIMHEADER_IMAGE = "community.wave.seqera.io/library/htslib_samtools_gawk:5366169842fe8538"

# Most of these are single-threaded I/O over modest inputs; sort is the one
# memory/throughput-sensitive step.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
SORT_RESOURCES = flyte.Resources(cpu=4, memory="12Gi")
# Upstream's process_low for the FASTQ conversions.
FASTQ_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
# calmd and markdup (upstream's process_medium), at the 4 CPUs upstream's tests run with.
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")


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


# As upstream: a BAM becomes a CRAM, anything else a BAM, then the output is
# indexed. The input and its index are linked side by side under their own
# names, as are the reference and its .fai (htslib may index the FASTA).
samtools_convert_cmd = shell.create(
    name="samtools_convert",
    image=SAMTOOLS_FASTQ_IMAGE,
    resources=FASTQ_RESOURCES,
    inputs={"alignment": File, "index": File | None, "fasta": File, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}}); FA=({{inputs.fasta}})
        shopt -s nullglob; IDX=({{inputs.index}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{IN[0]}}"); ln -s "${{IN[0]}}" "$NAME"
        if [ ${{#IDX[@]}} -gt 0 ]; then ln -s "${{IDX[0]}}" "$(basename "${{IDX[0]}}")"; fi
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
        if [ "${{NAME##*.}}" = bam ]; then EXT=cram; else EXT=bam; fi
        samtools view --threads $TASK_CPUS --reference "$FA_NAME" $ARGS "$NAME" -o "$OUT/$P.$EXT"
        samtools index -@$TASK_CPUS "$OUT/$P.$EXT"
    """,
)

# The FASTA is linked under its own name, which the dictionary is named after.
samtools_dict_cmd = shell.create(
    name="samtools_dict",
    image=SAMTOOLS_FASTQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$NAME"
        samtools dict $ARGS "$NAME" > {outputs.results}/"$NAME.dict"
    """,
)


# As upstream: `calmd` writes SAM (BAM with `-b` in args) to `<prefix>.bam`.
# The alignment, the FASTA and its .fai are linked side by side under their
# own names in a work dir, so the @PG command line names them as upstream's.
samtools_calmd_cmd = shell.create(
    name="samtools_calmd",
    image=SAMTOOLS_FASTQ_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"alignment": File, "fasta": File, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}}); FA=({{inputs.fasta}})
        shopt -s nullglob; FAI=({{inputs.fai}}); shopt -u nullglob
        ARGS={{inputs.args}}
        OUT={{outputs.results}}/{{inputs.prefix}}.bam
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{IN[0]}}"); ln -s "${{IN[0]}}" "$NAME"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
        samtools calmd -@ $TASK_CPUS $ARGS "$NAME" "$FA_NAME" > "$OUT"
    """,
)


def alignment_extension_script() -> str:
    """Bash setting ``EXT`` from ``--output-fmt`` in ``$ARGS`` (sam, bam, cram), else bam, as upstream."""
    return r"""
        case " $ARGS " in
            *"--output-fmt sam"*) EXT=sam ;;
            *"--output-fmt bam"*) EXT=bam ;;
            *"--output-fmt cram"*) EXT=cram ;;
            *) EXT=bam ;;
        esac
    """


# As upstream: the output format follows `--output-fmt` in args (default BAM),
# and `--threads` is the task's CPUs minus 1 (samtools counts extra threads).
# The reference, when given, is linked with its .fai into the work dir.
samtools_fixmate_cmd = shell.create(
    name="samtools_fixmate",
    image=SAMTOOLS_FASTQ_IMAGE,
    resources=FASTQ_RESOURCES,
    inputs={"alignment": File, "fasta": File | None, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        ARGS={{inputs.args}}
        {alignment_extension_script()}
        OUT={{outputs.results}}/{{inputs.prefix}}.$EXT
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{IN[0]}}"); ln -s "${{IN[0]}}" "$NAME"
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
            REF=(--reference "$FA_NAME")
        fi
        samtools fixmate $ARGS --threads $(( TASK_CPUS - 1 )) "${{REF[@]}}" "$NAME" "$OUT"
    """,
)

# As upstream: the output format follows `--output-fmt` in args (default BAM);
# temporary files go to `<prefix>.*` in the work dir (`-T`).
samtools_markdup_cmd = shell.create(
    name="samtools_markdup",
    image=SAMTOOLS_FASTQ_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"alignment": File, "fasta": File | None, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {alignment_extension_script()}
        OUT={{outputs.results}}/$P.$EXT
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{IN[0]}}"); ln -s "${{IN[0]}}" "$NAME"
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
            REF=(--reference "$FA_NAME")
        fi
        samtools markdup $ARGS "${{REF[@]}}" -@ $TASK_CPUS -T "$P" "$NAME" "$OUT"
    """,
)

# As upstream: the output format follows `--output-fmt` in args (default BAM),
# and the reads' flags follow their layout: two single-end files' first is
# passed bare and the rest with -0, a pair as -1/-2, one single-end file as
# -0, and one paired-end file as interleaved (-s). `read_names` ("/"-separated)
# keeps the caller's order, since a list[File] is globbed alphabetically.
samtools_import_cmd = shell.create(
    name="samtools_import",
    image=SAMTOOLS_FASTQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "read_names": str, "single_end": bool, "prefix": str, "args": str},
    defaults={"single_end": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R=({{inputs.reads}})
        NAMES={{inputs.read_names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {alignment_extension_script()}
        DIR=$(dirname "${{R[0]}}")
        IFS=/ read -ra ORDERED <<< "$NAMES"
        READS=()
        for NAME in "${{ORDERED[@]}}"; do READS+=("$DIR/$NAME"); done
        if [ ${{#READS[@]}} -gt 1 ] && [ {{inputs.single_end}} = true ]; then
            INPUT=("${{READS[0]}}")
            for f in "${{READS[@]:1}}"; do INPUT+=("-0$f"); done
        elif [ ${{#READS[@]}} -gt 1 ]; then
            INPUT=(-1 "${{READS[0]}}" -2 "${{READS[1]}}")
        elif [ {{inputs.single_end}} = true ]; then
            INPUT=(-0 "${{READS[0]}}")
        else
            INPUT=(-s "${{READS[0]}}")
        fi
        samtools import "${{INPUT[@]}}" $ARGS -@ $TASK_CPUS -o {{outputs.results}}/"$P.$EXT" || exit 1
    """,
)


# As upstream: idxstats (which needs the index beside the BAM) lists each
# reference's mapped and unmapped reads; the header loses the @SQ lines of
# references with neither, and the records are streamed back unchanged.
samtools_trimheader_cmd = shell.create(
    name="samtools_trimheader",
    image=SAMTOOLS_TRIMHEADER_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"bam": File, "bai": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        BAI=({{inputs.bai}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir input
        ln -s "${{BAM[0]}}" "input/$(basename "${{BAM[0]}}")"
        ln -s "${{BAI[0]}}" "input/$(basename "${{BAI[0]}}")"
        IN="input/$(basename "${{BAM[0]}}")"
        samtools idxstats --threads $TASK_CPUS "$IN" > idxstats.tsv || exit 1
        KEEP_USED='NR == FNR {{ if ($3 + $4 > 0) keep[$1]; next }} !/^@SQ/ || (substr($2, 4) in keep)'
        {{
            samtools view -H "$IN" | awk -F '\t' "$KEEP_USED" idxstats.tsv -
            samtools view --threads $TASK_CPUS "$IN"
        }} | samtools view $ARGS --threads $TASK_CPUS -b -o "$OUT/$P.bam" - || exit 1
        rm idxstats.tsv
    """,
)


# As upstream: split writes the mates, the reads flagged neither or both
# READ1 and READ2, and the singletons to their own gzipped files; otherwise
# the reads are written interleaved and bgzipped.
samtools_bam2fq_cmd = shell.create(
    name="samtools_bam2fq",
    image=SAMTOOLS_COLLATEFASTQ_IMAGE,
    resources=FASTQ_RESOURCES,
    inputs={"alignment": File, "prefix": str, "split": bool, "args": str},
    defaults={"split": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -o pipefail
        W=$(mktemp -d); cd "$W"
        if [ {{inputs.split}} = true ]; then
            samtools bam2fq $ARGS -@ $TASK_CPUS \
                -1 "$OUT/${{P}}_1.fq.gz" -2 "$OUT/${{P}}_2.fq.gz" \
                -0 "$OUT/${{P}}_other.fq.gz" -s "$OUT/${{P}}_singleton.fq.gz" \
                "${{IN[0]}}" || exit 1
        else
            samtools bam2fq $ARGS -@ $TASK_CPUS "${{IN[0]}}" | bgzip > "$OUT/${{P}}_interleaved.fq.gz" || exit 1
        fi
    """,
)

# As upstream: the reference and the intervals (`-l`) are optional; the
# alignment is linked beside its index (needed for a region, `-r` in args),
# and the pileup is bgzipped.
samtools_mpileup_cmd = shell.create(
    name="samtools_mpileup",
    image=SAMTOOLS_COLLATEFASTQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "alignment": File,
        "index": File | None,
        "intervals": File | None,
        "fasta": File | None,
        "fai": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.alignment}})
        shopt -s nullglob
        IDX=({{inputs.index}}); BED=({{inputs.intervals}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{IN[0]}}"); ln -s "${{IN[0]}}" "$NAME"
        if [ ${{#IDX[@]}} -gt 0 ]; then ln -s "${{IDX[0]}}" "$(basename "${{IDX[0]}}")"; fi
        OPTIONAL=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_LINK=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_LINK"
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_LINK.fai"; fi
            OPTIONAL+=(--fasta-ref "$FA_LINK")
        fi
        INTERVALS=()
        if [ ${{#BED[@]}} -gt 0 ]; then INTERVALS=(-l "${{BED[0]}}"); fi
        samtools mpileup "${{OPTIONAL[@]}}" --output "$P.mpileup" $ARGS "${{INTERVALS[@]}}" "$NAME" || exit 1
        bgzip "$P.mpileup" || exit 1
        mv "$P.mpileup.gz" {{outputs.results}}/
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
samtools_fastq_env = flyte.TaskEnvironment.from_task(
    "samtools_fastq",
    samtools_fastq_cmd.as_task(),
    samtools_convert_cmd.as_task(),
    samtools_dict_cmd.as_task(),
    samtools_calmd_cmd.as_task(),
    samtools_fixmate_cmd.as_task(),
    samtools_markdup_cmd.as_task(),
    samtools_import_cmd.as_task(),
)
samtools_collatefastq_env = flyte.TaskEnvironment.from_task(
    "samtools_collatefastq",
    samtools_collatefastq_cmd.as_task(),
    samtools_bam2fq_cmd.as_task(),
    samtools_mpileup_cmd.as_task(),
)
samtools_trimheader_env = flyte.TaskEnvironment.from_task("samtools_trimheader", samtools_trimheader_cmd.as_task())
env = flyte.TaskEnvironment(
    name="samtools",
    depends_on=[samtools_tools_env, samtools_fastq_env, samtools_collatefastq_env, samtools_trimheader_env],
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


@dataclass
class SamtoolsConvertResult:
    alignment: File  # <prefix>.cram from a BAM, else <prefix>.bam
    index: File  # <prefix>.cram.crai or <prefix>.bam.bai
    results: Dir  # everything samtools wrote


async def samtools_convert(
    alignment: File, fasta: File, prefix: str, index: File | None = None, fai: File | None = None, args: str = ""
) -> SamtoolsConvertResult:
    """Convert a BAM ``alignment`` to ``<prefix>.cram`` (anything else to ``<prefix>.bam``) and index it.

    ``fasta`` (with its ``fai``) is the reference, which a CRAM needs.
    """
    results = await samtools_convert_cmd(
        alignment=alignment, index=index, fasta=fasta, fai=fai, prefix=prefix, args=args
    )
    extension, index_suffix = ("cram", "crai") if alignment.name.endswith(".bam") else ("bam", "bai")
    converted = await results.get_file(f"{prefix}.{extension}")
    converted_index = await results.get_file(f"{prefix}.{extension}.{index_suffix}")
    if converted is None or converted_index is None:
        raise FileNotFoundError(f"samtools convert wrote no {prefix}.{extension} and its .{index_suffix}")
    return SamtoolsConvertResult(alignment=converted, index=converted_index, results=results)


async def samtools_dict(fasta: File, args: str = "") -> File:
    """``<fasta name>.dict``: the sequence dictionary of ``fasta``."""
    results = await samtools_dict_cmd(fasta=fasta, args=args)
    dictionary = await results.get_file(f"{fasta.name}.dict")
    if dictionary is None:
        raise FileNotFoundError(f"samtools dict wrote no {fasta.name}.dict")
    return dictionary


async def pick_alignment(results: Dir, prefix: str, args: str, tool: str) -> File:
    """``<prefix>.<sam|bam|cram>`` from ``results``, the extension following ``--output-fmt`` in ``args``."""
    extension = "bam"
    for fmt in ("sam", "bam", "cram"):
        if f"--output-fmt {fmt}" in args:
            extension = fmt
            break
    found = await results.get_file(f"{prefix}.{extension}")
    if found is None:
        raise FileNotFoundError(f"samtools {tool} wrote no {prefix}.{extension}")
    return found


async def samtools_calmd(alignment: File, fasta: File, prefix: str, fai: File | None = None, args: str = "") -> File:
    """``alignment`` with its MD and NM tags recalculated against ``fasta``, as ``<prefix>.bam``.

    As upstream, the output is SAM text unless ``args`` asks for BAM (``-b``).
    """
    if alignment.name == f"{prefix}.bam":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await samtools_calmd_cmd(alignment=alignment, fasta=fasta, fai=fai, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.bam")
    if out is None:
        raise FileNotFoundError(f"samtools calmd wrote no {prefix}.bam")
    return out


async def samtools_fixmate(
    alignment: File, prefix: str, fasta: File | None = None, fai: File | None = None, args: str = ""
) -> File:
    """A name-grouped ``alignment`` with its mate fields filled in, as ``<prefix>.<ext>``.

    The extension follows ``--output-fmt`` in ``args`` (default bam); a CRAM needs ``fasta``.
    Upstream's default prefix is ``<id>_fixmate``.
    """
    if alignment.name in (f"{prefix}.sam", f"{prefix}.bam", f"{prefix}.cram"):
        raise ValueError("input and output names are the same; choose another prefix")
    results = await samtools_fixmate_cmd(alignment=alignment, fasta=fasta, fai=fai, prefix=prefix, args=args)
    return await pick_alignment(results, prefix, args, "fixmate")


async def samtools_markdup(
    alignment: File, prefix: str, fasta: File | None = None, fai: File | None = None, args: str = ""
) -> File:
    """A coordinate-sorted ``alignment`` with its duplicates marked, as ``<prefix>.<ext>``.

    The extension follows ``--output-fmt`` in ``args`` (default bam); a CRAM needs ``fasta``.
    """
    if alignment.name in (f"{prefix}.sam", f"{prefix}.bam", f"{prefix}.cram"):
        raise ValueError("input and output names are the same; choose another prefix")
    results = await samtools_markdup_cmd(alignment=alignment, fasta=fasta, fai=fai, prefix=prefix, args=args)
    return await pick_alignment(results, prefix, args, "markdup")


async def samtools_import(reads: list[File], prefix: str, single_end: bool = False, args: str = "") -> File:
    """FASTQ ``reads`` as an unaligned ``<prefix>.<ext>`` (``--output-fmt`` in ``args``, default bam).

    ``reads`` is one single-end file, several single-end files, a pair (R1, R2) or, with
    ``single_end`` false, one interleaved file. Their names must be distinct.
    """
    if not reads:
        raise ValueError("samtools import needs at least one FASTQ")
    names = [f.name for f in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"samtools import needs distinct FASTQ names, got {names}")
    results = await samtools_import_cmd(
        reads=reads, read_names="/".join(names), single_end=single_end, prefix=prefix, args=args
    )
    return await pick_alignment(results, prefix, args, "import")


async def samtools_trimheader(bam: File, bai: File, prefix: str, args: str = "") -> File:
    """``bam`` (indexed by ``bai``) with only the @SQ lines of references with reads, as ``<prefix>.bam``."""
    if bam.name == f"{prefix}.bam":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await samtools_trimheader_cmd(bam=bam, bai=bai, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.bam")
    if out is None:
        raise FileNotFoundError(f"samtools trimheader wrote no {prefix}.bam")
    return out


@dataclass
class SamtoolsBam2fqResult:
    fastq: list[File]  # <prefix>_1.fq.gz and <prefix>_2.fq.gz when split, else <prefix>_interleaved.fq.gz
    other: File | None  # <prefix>_other.fq.gz, when split
    singleton: File | None  # <prefix>_singleton.fq.gz, when split


async def samtools_bam2fq(alignment: File, prefix: str, split: bool = False, args: str = "") -> SamtoolsBam2fqResult:
    """The reads of ``alignment`` as gzipped FASTQ: split by mate, or interleaved."""
    results = await samtools_bam2fq_cmd(alignment=alignment, prefix=prefix, split=split, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"samtools bam2fq wrote no {name}")
        return found

    if not split:
        return SamtoolsBam2fqResult(fastq=[await pick(f"{prefix}_interleaved.fq.gz")], other=None, singleton=None)
    return SamtoolsBam2fqResult(
        fastq=[await pick(f"{prefix}_1.fq.gz"), await pick(f"{prefix}_2.fq.gz")],
        other=await pick(f"{prefix}_other.fq.gz"),
        singleton=await pick(f"{prefix}_singleton.fq.gz"),
    )


async def samtools_mpileup(
    alignment: File,
    prefix: str,
    index: File | None = None,
    intervals: File | None = None,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> File:
    """The pileup of ``alignment``: ``<prefix>.mpileup.gz``.

    ``index`` is needed for a region (``-r`` in ``args``); ``intervals`` is a
    BED of positions (``-l``).
    """
    results = await samtools_mpileup_cmd(
        alignment=alignment, index=index, intervals=intervals, fasta=fasta, fai=fai, prefix=prefix, args=args
    )
    found = await results.get_file(f"{prefix}.mpileup.gz")
    if found is None:
        raise FileNotFoundError(f"samtools mpileup wrote no {prefix}.mpileup.gz")
    return found
