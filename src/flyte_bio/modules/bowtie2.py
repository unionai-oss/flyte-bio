"""bowtie2 — Bowtie 2 short-read alignment.

Exposes:

- :func:`bowtie2_build` — build a Bowtie 2 index of a FASTA. The index files
  (``<name>.{1,2,3,4}.bt2`` and ``<name>.rev.{1,2}.bt2``, or ``.bt2l`` with
  ``--large-index``, named after the FASTA without its extension, as upstream)
  land in one output Dir.
- :func:`bowtie2_align` — align one sample's reads (one or two mates), piped
  into ``samtools sort`` (or ``view`` with ``sort_bam=False``).

Outputs of :func:`bowtie2_align` keep upstream's names:

- the alignment: ``<prefix>.bam``, or ``<prefix>.sam`` / ``<prefix>.cram`` when
  ``args2`` sets ``--output-fmt`` / ``-O``; a CRAM needs ``fasta``;
- its index, when ``args2`` has ``--write-index``: ``<prefix>.bam.csi`` or
  ``<prefix>.cram.crai``;
- Bowtie 2's alignment summary: ``<prefix>.bowtie2.log``;
- with ``save_unaligned``: ``<prefix>.unmapped.fastq.gz`` (single-end), or the
  pairs that didn't align concordantly, ``<prefix>.unmapped_1.fastq.gz`` and
  ``<prefix>.unmapped_2.fastq.gz``.

As upstream, reads get the read group ``--rg-id <prefix> --rg SM:<prefix>``
unless ``args`` sets ``--rg-id``, and Bowtie 2 and samtools both use the task's
CPUs.
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BOWTIE2_IMAGE = "community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:dbd18842ad7cb278"
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")

# The output format, as upstream reads it from samtools' arguments.
OUTPUT_FORMAT = r"(--output-fmt|-O)+\s+(\S+)"

bowtie2_build_cmd = shell.create(
    name="bowtie2_build",
    image=BOWTIE2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        NAME=$(basename "${{FA[0]}}")
        ARGS={{inputs.args}}
        bowtie2-build $ARGS --threads $TASK_CPUS "${{FA[0]}}" {{outputs.index}}/"${{NAME%.*}}"
    """,
)

# As upstream: the index is found by its .rev.1.bt2 (or .bt2l) file, and
# Bowtie 2's summary goes to <prefix>.bowtie2.log. Bowtie 2 names the
# unaligned pairs <prefix>.unmapped.fastq.{1,2}.gz; they're renamed to
# <prefix>.unmapped_{1,2}.fastq.gz.
bowtie2_align_cmd = shell.create(
    name="bowtie2_align",
    image=BOWTIE2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "fasta": File | None,
        "prefix": str,
        "save_unaligned": bool,
        "sort_bam": bool,
        "args": str,
        "args2": str,
    },
    defaults={"save_unaligned": False, "sort_bam": True, "args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); FA=({{inputs.fasta}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}

        INDEX=$(find -L {{inputs.index}} -name "*.rev.1.bt2" | sed 's/\.rev\.1\.bt2$//')
        if [ -z "$INDEX" ]; then
            INDEX=$(find -L {{inputs.index}} -name "*.rev.1.bt2l" | sed 's/\.rev\.1\.bt2l$//')
        fi
        if [ -z "$INDEX" ]; then echo "Bowtie2 index files not found" >&2; exit 1; fi

        FORMAT='(--output-fmt|-O)+[[:space:]]+([^[:space:]]+)'
        EXT=bam
        if [[ $ARGS2 =~ $FORMAT ]]; then EXT=${{BASH_REMATCH[2],,}}; fi
        REF=()
        if [ $EXT = cram ]; then
            if [ ${{#FA[@]}} -eq 0 ]; then echo "Fasta reference is required for CRAM output" >&2; exit 1; fi
            W=$(mktemp -d)
            ln -s "${{FA[0]}}" "$W/$(basename "${{FA[0]}}")"
            REF=(--reference "$W/$(basename "${{FA[0]}}")")
        fi

        RG=(--rg-id "$P" --rg "SM:$P")
        if [[ "$ARGS" == *--rg-id* ]]; then RG=(); fi
        if [ ${{#R2[@]}} -eq 0 ]; then
            READS=(-U "${{R1[0]}}")
            UNALIGNED=(--un-gz "$OUT/$P.unmapped.fastq.gz")
        else
            READS=(-1 "${{R1[0]}}" -2 "${{R2[0]}}")
            UNALIGNED=(--un-conc-gz "$OUT/$P.unmapped.fastq.gz")
        fi
        if [ {{inputs.save_unaligned}} != true ]; then UNALIGNED=(); fi
        SAMTOOLS=view
        if [ {{inputs.sort_bam}} = true ]; then SAMTOOLS=sort; fi

        bowtie2 -x "$INDEX" "${{READS[@]}}" --threads $TASK_CPUS "${{UNALIGNED[@]}}" "${{RG[@]}}" $ARGS \
            2> "$OUT/$P.bowtie2.log" \
            | samtools $SAMTOOLS $ARGS2 --threads $TASK_CPUS "${{REF[@]}}" -o "$OUT/$P.$EXT" -
        cat "$OUT/$P.bowtie2.log" >&2

        for MATE in 1 2; do
            if [ -f "$OUT/$P.unmapped.fastq.$MATE.gz" ]; then
                mv "$OUT/$P.unmapped.fastq.$MATE.gz" "$OUT/${{P}}.unmapped_$MATE.fastq.gz"
            fi
        done
        true
    """,
)

env = flyte.TaskEnvironment.from_task(
    "bowtie2",
    bowtie2_build_cmd.as_task(),
    bowtie2_align_cmd.as_task(),
)


async def bowtie2_build(fasta: File, args: str = "") -> Dir:
    """Bowtie 2 index of ``fasta``: a Dir of ``<fasta name>.{1,2,3,4,rev.1,rev.2}.bt2``."""
    return await bowtie2_build_cmd(fasta=fasta, args=args)


@dataclass
class Bowtie2AlignResult:
    alignment: File  # <prefix>.bam, .sam or .cram
    alignment_index: File | None  # <prefix>.bam.csi or <prefix>.cram.crai, with --write-index in args2
    log: File  # <prefix>.bowtie2.log
    unaligned: list[File]  # with save_unaligned: one file (single-end) or mates 1 and 2
    results: Dir  # everything the task wrote


async def bowtie2_align(
    reads_1: File,
    reads_2: File | None,
    index: Dir,
    prefix: str,
    fasta: File | None = None,
    save_unaligned: bool = False,
    sort_bam: bool = True,
    args: str = "",
    args2: str = "",
) -> Bowtie2AlignResult:
    """Align one sample's reads with Bowtie 2; ``args`` go to Bowtie 2, ``args2`` to samtools.

    ``reads_2`` is None for single-end reads. ``fasta`` is the reference, needed
    for CRAM output (``--output-fmt cram`` in ``args2``).
    """
    found = re.search(OUTPUT_FORMAT, args2)
    extension = found.group(2).lower() if found else "bam"
    if extension == "cram" and fasta is None:
        raise ValueError("Fasta reference is required for CRAM output")
    results = await bowtie2_align_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        fasta=fasta,
        prefix=prefix,
        save_unaligned=save_unaligned,
        sort_bam=sort_bam,
        args=args,
        args2=args2,
    )

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"bowtie2 align wrote no {name}")
        return file

    alignment_index = None
    for suffix in ("csi", "crai"):
        alignment_index = alignment_index or await results.get_file(f"{prefix}.{extension}.{suffix}")
    unaligned = []
    if save_unaligned:
        mates = ["_1", "_2"] if reads_2 is not None else [""]
        unaligned = [await pick(f"{prefix}.unmapped{mate}.fastq.gz") for mate in mates]
    return Bowtie2AlignResult(
        alignment=await pick(f"{prefix}.{extension}"),
        alignment_index=alignment_index,
        log=await pick(f"{prefix}.bowtie2.log"),
        unaligned=unaligned,
        results=results,
    )
