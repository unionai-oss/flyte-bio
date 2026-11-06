"""minimap2 — alignment of long or short reads, and of assemblies, to a reference.

Exposes:

- :func:`minimap2_index` — build a minimap2 index of a FASTA:
  ``<fasta name without its last extension>.mmi``, as upstream.
- :func:`minimap2_align` — align one sample's reads to a reference FASTA or
  ``.mmi`` index. With ``bam_format`` the alignments are piped into
  ``samtools sort`` to ``<prefix>.bam`` (indexed to ``<prefix>.bam.<ext>``
  with ``bam_index_extension``, e.g. ``bai`` or ``csi``); otherwise minimap2
  writes ``<prefix>.paf`` (with CIGARs, ``-c``, when ``cigar_paf_format``).

As upstream:

- ``reads`` is one or two FASTQs, or one SAM/BAM/CRAM. Alignments are first
  turned back into reads (``samtools reset | samtools fastq``) and need a
  reference.
- Without a reference the reads are aligned to themselves.
- ``cigar_bam`` adds ``-L``, which moves CIGARs with more than 65535
  operations into the ``CG`` tag of BAM output.
- minimap2 gets the task's CPUs as threads; samtools gets one fewer.
- Inputs are linked into a working directory under their own names, so the
  BAM header's ``@PG`` command lines name them as upstream's do.

``args`` go to minimap2, ``args2`` to ``samtools sort``, ``args3`` to
``samtools reset`` and ``args4`` to ``samtools fastq``.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MINIMAP2_IMAGE = "community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd"
# Alignment adds samtools, for sorting and for alignment input.
MINIMAP2_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"

# Upstream labels indexing process_low and alignment process_high. A human
# index takes ~11 GB of memory to build or load.
INDEX_RESOURCES = flyte.Resources(cpu=2, memory="16Gi")
ALIGN_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

ALIGNMENT_EXTENSIONS = ("sam", "bam", "cram")

minimap2_index_cmd = shell.create(
    name="minimap2_index",
    image=MINIMAP2_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        FA=({{inputs.fasta}})
        NAME=$(basename "${{FA[0]}}")
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{FA[0]}}" "$NAME"
        minimap2 -t $TASK_CPUS -d "${{NAME%.*}}.mmi" $ARGS "$NAME"
        mv ./*.mmi {{outputs.results}}/
    """,
)

minimap2_align_cmd = shell.create(
    name="minimap2_align",
    image=MINIMAP2_SAMTOOLS_IMAGE,
    resources=ALIGN_RESOURCES,
    inputs={
        "reads": list[File],
        "reference": File | None,
        "prefix": str,
        "bam_format": bool,
        "bam_index_extension": str,
        "cigar_paf_format": bool,
        "cigar_bam": bool,
        "args": str,
        "args2": str,
        "args3": str,
        "args4": str,
    },
    defaults={
        "bam_format": True,
        "bam_index_extension": "",
        "cigar_paf_format": False,
        "cigar_bam": False,
        "args": "",
        "args2": "",
        "args3": "",
        "args4": "",
    },
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        READS=({{inputs.reads}})
        shopt -s nullglob; REF=({{inputs.reference}}); shopt -u nullglob
        P={{inputs.prefix}}
        INDEX_EXT={{inputs.bam_index_extension}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        ARGS3={{inputs.args3}}
        ARGS4={{inputs.args4}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        NAMES=()
        for f in "${{READS[@]}}"; do ln -s "$f" "$(basename "$f")"; NAMES+=("$(basename "$f")"); done
        ALIGNMENT_INPUT=false
        case "${{NAMES[0]##*.}}" in sam|bam|cram) ALIGNMENT_INPUT=true ;; esac
        TARGET=("${{NAMES[@]}}")
        if [ ${{#REF[@]}} -gt 0 ]; then
            ln -s "${{REF[0]}}" "$(basename "${{REF[0]}}")"
            TARGET=("$(basename "${{REF[0]}}")")
        elif [ $ALIGNMENT_INPUT = true ]; then
            echo "minimap2_align: alignment input needs a reference" >&2
            exit 1
        fi
        QUERY=("${{NAMES[@]}}")
        if [ $ALIGNMENT_INPUT = true ]; then QUERY=(-); fi
        SAMTOOLS_THREADS=$(( TASK_CPUS - 1 ))

        BAM_FORMAT={{inputs.bam_format}}
        CIGAR=()
        if [ $BAM_FORMAT = true ]; then
            if [ {{inputs.cigar_bam}} = true ]; then CIGAR=(-L); fi
        elif [ {{inputs.cigar_paf_format}} = true ]; then
            CIGAR=(-c)
        fi
        SORT_OUT=(-o "$P.bam")
        if [ -n "$INDEX_EXT" ]; then SORT_OUT=(-o "$P.bam##idx##$P.bam.$INDEX_EXT" --write-index); fi

        align() {{
            if [ $BAM_FORMAT = true ]; then
                minimap2 $ARGS -t $TASK_CPUS "${{TARGET[@]}}" "${{QUERY[@]}}" "${{CIGAR[@]}}" -a \
                    | samtools sort -@ $SAMTOOLS_THREADS "${{SORT_OUT[@]}}" $ARGS2
            else
                minimap2 $ARGS -t $TASK_CPUS "${{TARGET[@]}}" "${{QUERY[@]}}" "${{CIGAR[@]}}" -o "$P.paf"
            fi
        }}
        if [ $ALIGNMENT_INPUT = true ]; then
            samtools reset --threads $SAMTOOLS_THREADS $ARGS3 "${{NAMES[0]}}" \
                | samtools fastq --threads $SAMTOOLS_THREADS $ARGS4 \
                | align
        else
            align
        fi
        for f in "$P.paf" "$P.bam" "$P.bam.$INDEX_EXT"; do
            if [ -f "$f" ]; then mv "$f" "$OUT"/; fi
        done
    """,
)

# One image per environment: alignment's adds samtools, so `env` aggregates two.
minimap2_env = flyte.TaskEnvironment.from_task("minimap2_tools", minimap2_index_cmd.as_task())
minimap2_samtools_env = flyte.TaskEnvironment.from_task("minimap2_samtools", minimap2_align_cmd.as_task())
env = flyte.TaskEnvironment(name="minimap2", depends_on=[minimap2_env, minimap2_samtools_env])


@dataclass
class Minimap2AlignResult:
    bam: File | None  # <prefix>.bam, with bam_format
    index: File | None  # <prefix>.bam.<bam_index_extension>, with bam_format and an extension
    paf: File | None  # <prefix>.paf, without bam_format
    results: Dir  # everything minimap2 wrote


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def minimap2_index(fasta: File, args: str = "") -> File:
    """minimap2 index of ``fasta``: ``<fasta name without its last extension>.mmi``."""
    results = await minimap2_index_cmd(fasta=fasta, args=args)
    name = file_name(fasta)
    stem = name.rsplit(".", 1)[0] if "." in name else name
    index = await results.get_file(f"{stem}.mmi")
    if index is None:
        raise FileNotFoundError(f"minimap2 index wrote no {stem}.mmi")
    return index


async def minimap2_align(
    reads: list[File],
    reference: File | None,
    prefix: str = "sample",
    bam_format: bool = True,
    bam_index_extension: str = "",
    cigar_paf_format: bool = False,
    cigar_bam: bool = False,
    args: str = "",
    args2: str = "",
    args3: str = "",
    args4: str = "",
) -> Minimap2AlignResult:
    """Align ``reads`` (one or two FASTQs, or one SAM/BAM/CRAM) to ``reference`` (a FASTA or ``.mmi``).

    Without a reference, FASTQ reads are aligned to themselves; alignment
    input needs one.
    """
    if not reads:
        raise ValueError("minimap2_align needs reads")
    if reference is None and file_name(reads[0]).rsplit(".", 1)[-1] in ALIGNMENT_EXTENSIONS:
        raise ValueError("minimap2_align: alignment input needs a reference")
    results = await minimap2_align_cmd(
        reads=reads,
        reference=reference,
        prefix=prefix,
        bam_format=bam_format,
        bam_index_extension=bam_index_extension,
        cigar_paf_format=cigar_paf_format,
        cigar_bam=cigar_bam,
        args=args,
        args2=args2,
        args3=args3,
        args4=args4,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"minimap2 align wrote no {name}")
        return found

    return Minimap2AlignResult(
        bam=await pick(f"{prefix}.bam") if bam_format else None,
        index=await pick(f"{prefix}.bam.{bam_index_extension}") if bam_format and bam_index_extension else None,
        paf=None if bam_format else await pick(f"{prefix}.paf"),
        results=results,
    )
