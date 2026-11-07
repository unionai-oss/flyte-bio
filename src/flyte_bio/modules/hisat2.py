"""hisat2 — HISAT2 splice-aware short-read alignment.

Exposes:

- :func:`hisat2_extractsplicesites` — the known splice sites of a GTF,
  ``<gtf name>.splice_sites.txt`` (HISAT2's ``hisat2_extract_splice_sites.py``).
- :func:`hisat2_build` — build a HISAT2 index of a FASTA. The index files
  (``<name>.{1..8}.ht2``, named after the FASTA without its extension, as
  upstream) land in one output Dir. With a GTF and splice sites, the index is
  splice-aware, but only when the task has at least ``splice_aware_memory_gb``
  GiB of memory, as upstream decides from ``task.memory``.
- :func:`hisat2_align` — align one sample's reads (one or two mates), piped
  into ``samtools view``, keeping only mapped primary alignments (and, for
  pairs, only concordant pairs with both mates mapped).

Outputs of :func:`hisat2_align` keep upstream's names:

- the alignment, unsorted: ``<prefix>.bam``;
- HISAT2's alignment summary: ``<prefix>.hisat2.summary.log``;
- with ``save_unaligned``: ``<prefix>.unmapped.fastq.gz`` (single-end), or the
  pairs that didn't align concordantly, ``<prefix>.unmapped_1.fastq.gz`` and
  ``<prefix>.unmapped_2.fastq.gz``.

As upstream, reads get the read group ``--rg-id <prefix> --rg SM:<prefix>``
unless ``args`` sets ``--rg-id``, and HISAT2 uses the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream builds with the plain HISAT2 image and aligns with one that adds samtools.
HISAT2_IMAGE = "community.wave.seqera.io/library/hisat2:2.2.3--2616fa83d3b9d8f8"
HISAT2_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/hisat2_samtools:a0c9b8ccf8116a89"
# A splice-aware index of a human genome needs far more memory (upstream
# pipelines default to a 200 GB threshold for it); override the resources then.
BUILD_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")
ALIGN_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")
EXTRACT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

hisat2_extractsplicesites_cmd = shell.create(
    name="hisat2_extractsplicesites",
    image=HISAT2_SAMTOOLS_IMAGE,
    resources=EXTRACT_RESOURCES,
    inputs={"gtf": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        NAME=$(basename "${GTF[0]}")
        ARGS={inputs.args}
        hisat2_extract_splice_sites.py $ARGS "${GTF[0]}" > {outputs.results}/"${NAME%.*}".splice_sites.txt || exit 1
    """,
)

# As upstream: splice sites and exons only go into the index when the task's
# memory (whole GiB) reaches the threshold; 0 never uses them.
hisat2_build_cmd = shell.create(
    name="hisat2_build",
    image=HISAT2_IMAGE,
    resources=BUILD_RESOURCES,
    inputs={
        "fasta": File,
        "gtf": File | None,
        "splicesites": File | None,
        "splice_aware_memory_gb": int,
        "args": str,
    },
    defaults={"splice_aware_memory_gb": 0, "args": ""},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; GTF=({{inputs.gtf}}); SS=({{inputs.splicesites}}); shopt -u nullglob
        NAME=$(basename "${{FA[0]}}")
        ARGS={{inputs.args}}
        THRESHOLD={{inputs.splice_aware_memory_gb}}
        AVAILABLE=$(( TASK_MEMORY_MB / 1024 ))

        SPLICE=()
        if [ "$THRESHOLD" -gt 0 ] && [ "$AVAILABLE" -ge "$THRESHOLD" ]; then
            echo "[HISAT2 index build] $AVAILABLE GB available, using splice sites and exons" >&2
            if [ ${{#SS[@]}} -gt 0 ]; then SPLICE+=(--ss "${{SS[0]}}"); fi
            if [ ${{#GTF[@]}} -gt 0 ]; then
                W=$(mktemp -d)
                GTF_NAME=$(basename "${{GTF[0]}}")
                EXONS="$W/${{GTF_NAME%.*}}.exons.txt"
                hisat2_extract_exons.py "${{GTF[0]}}" > "$EXONS" || exit 1
                SPLICE+=(--exon "$EXONS")
            fi
        else
            echo "[HISAT2 index build] Only $AVAILABLE GB available (< $THRESHOLD GB), so NOT splice-aware" >&2
        fi
        hisat2-build -p $TASK_CPUS "${{SPLICE[@]}}" $ARGS "${{FA[0]}}" {{outputs.index}}/"${{NAME%.*}}" || exit 1
    """,
)

# As upstream: the index is found by its .1.ht2 (or .1.ht2l) file. HISAT2
# names the unaligned pairs <prefix>.unmapped.fastq.{1,2}.gz; they're renamed
# to <prefix>.unmapped_{1,2}.fastq.gz. HISAT2's wrapper is Perl, which exits
# 255 on errors; any failure exits 1 instead.
hisat2_align_cmd = shell.create(
    name="hisat2_align",
    image=HISAT2_SAMTOOLS_IMAGE,
    resources=ALIGN_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "splicesites": File | None,
        "prefix": str,
        "save_unaligned": bool,
        "args": str,
    },
    defaults={"save_unaligned": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); SS=({{inputs.splicesites}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}

        INDEX=$(find -L {{inputs.index}} -name "*.1.ht2*" | sed 's/\.1\.ht2.*$//')
        if [ -z "$INDEX" ]; then echo "HISAT2 index files not found" >&2; exit 1; fi

        SPLICE=()
        if [ ${{#SS[@]}} -gt 0 ]; then SPLICE=(--known-splicesite-infile "${{SS[0]}}"); fi
        RG=(--rg-id "$P" --rg "SM:$P")
        if [[ "$ARGS" == *--rg-id* ]]; then RG=(); fi
        if [ ${{#R2[@]}} -eq 0 ]; then
            READS=(-U "${{R1[0]}}")
            UNALIGNED=(--un-gz "$OUT/$P.unmapped.fastq.gz")
            PAIRED=()
            FILTER=(-F 4 -F 256)
        else
            READS=(-1 "${{R1[0]}}" -2 "${{R2[0]}}")
            UNALIGNED=(--un-conc-gz "$OUT/$P.unmapped.fastq.gz")
            PAIRED=(--no-mixed --no-discordant)
            FILTER=(-F 4 -F 8 -F 256)
        fi
        if [ {{inputs.save_unaligned}} != true ]; then UNALIGNED=(); fi

        hisat2 -x "$INDEX" "${{READS[@]}}" "${{SPLICE[@]}}" --summary-file "$OUT/$P.hisat2.summary.log" \
            --threads $TASK_CPUS "${{RG[@]}}" "${{UNALIGNED[@]}}" "${{PAIRED[@]}}" $ARGS \
            | samtools view -bS "${{FILTER[@]}}" - > "$OUT/$P.bam" || exit 1

        for MATE in 1 2; do
            if [ -f "$OUT/$P.unmapped.fastq.$MATE.gz" ]; then
                mv "$OUT/$P.unmapped.fastq.$MATE.gz" "$OUT/${{P}}.unmapped_$MATE.fastq.gz"
            fi
        done
        true
    """,
)

hisat2_build_env = flyte.TaskEnvironment.from_task("hisat2_build", hisat2_build_cmd.as_task())
hisat2_samtools_env = flyte.TaskEnvironment.from_task(
    "hisat2_samtools",
    hisat2_extractsplicesites_cmd.as_task(),
    hisat2_align_cmd.as_task(),
)
env = flyte.TaskEnvironment(name="hisat2", depends_on=[hisat2_build_env, hisat2_samtools_env])


async def hisat2_extractsplicesites(gtf: File, args: str = "") -> File:
    """The known splice sites of ``gtf``: ``<gtf name without extension>.splice_sites.txt``."""
    results = await hisat2_extractsplicesites_cmd(gtf=gtf, args=args)
    name = gtf.name.rsplit(".", 1)[0] + ".splice_sites.txt"
    file = await results.get_file(name)
    if file is None:
        raise FileNotFoundError(f"hisat2_extract_splice_sites.py wrote no {name}")
    return file


async def hisat2_build(
    fasta: File,
    gtf: File | None = None,
    splicesites: File | None = None,
    splice_aware_memory_gb: int | None = None,
    args: str = "",
) -> Dir:
    """HISAT2 index of ``fasta``: a Dir of ``<fasta name>.{1..8}.ht2``.

    ``gtf`` (its exons) and ``splicesites`` go into the index only when the
    task has at least ``splice_aware_memory_gb`` GiB of memory; with None, as
    upstream's ``false``, never.
    """
    return await hisat2_build_cmd(
        fasta=fasta,
        gtf=gtf,
        splicesites=splicesites,
        splice_aware_memory_gb=splice_aware_memory_gb or 0,
        args=args,
    )


@dataclass
class Hisat2AlignResult:
    bam: File  # <prefix>.bam
    summary: File  # <prefix>.hisat2.summary.log
    unaligned: list[File]  # with save_unaligned: one file (single-end) or mates 1 and 2
    results: Dir  # everything the task wrote


async def hisat2_align(
    reads_1: File,
    reads_2: File | None,
    index: Dir,
    prefix: str,
    splicesites: File | None = None,
    save_unaligned: bool = False,
    args: str = "",
) -> Hisat2AlignResult:
    """Align one sample's reads with HISAT2; ``reads_2`` is None for single-end reads.

    ``splicesites`` (from :func:`hisat2_extractsplicesites`) are passed as known
    splice sites.
    """
    results = await hisat2_align_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        splicesites=splicesites,
        prefix=prefix,
        save_unaligned=save_unaligned,
        args=args,
    )

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"hisat2 align wrote no {name}")
        return file

    unaligned = []
    if save_unaligned:
        mates = ["_1", "_2"] if reads_2 is not None else [""]
        unaligned = [await pick(f"{prefix}.unmapped{mate}.fastq.gz") for mate in mates]
    return Hisat2AlignResult(
        bam=await pick(f"{prefix}.bam"),
        summary=await pick(f"{prefix}.hisat2.summary.log"),
        unaligned=unaligned,
        results=results,
    )
