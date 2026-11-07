"""dragmap — DRAGMAP, Illumina's open-source DRAGEN short-read mapper.

Exposes:

- :func:`dragmap_hashtable` — build DRAGMAP's hash table of a FASTA. Its files
  (``hash_table.cmp``, ``hash_table.cfg``, ``reference.bin`` and others) land in
  one output Dir.
- :func:`dragmap_align` — align one sample's reads (one or two mates), piped
  into ``samtools sort`` (or ``view`` with ``sort_bam=False``).

Outputs of :func:`dragmap_align` keep upstream's names:

- the alignment: ``<prefix>.bam``, or ``<prefix>.sam`` / ``<prefix>.cram`` when
  ``args2`` sets ``--output-fmt`` / ``-O``; a CRAM needs ``fasta``;
- its index, when ``args2`` has ``--write-index``: ``<prefix>.bam.csi`` or
  ``<prefix>.cram.crai``;
- DRAGMAP's log (its stderr): ``<prefix>.dragmap.log``.

As upstream, DRAGMAP and samtools both use the task's CPUs. Upstream pins
DRAGMAP 1.2.1 (1.3.0 has an open alignment bug).
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DRAGMAP_IMAGE = "quay.io/biocontainers/dragmap:1.2.1--h72d16da_1"
# DRAGMAP with samtools and pigz.
DRAGMAP_SAMTOOLS_IMAGE = (
    "quay.io/biocontainers/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d"
    ":df80ed8d23d0a2c43181a2b3dd1b39f2d00fab5c-0"
)
# The hash table of a human genome takes far more memory to build and to align
# against (upstream's process_high gives 72 GB); override the resources then.
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")

# The output format, as upstream reads it from samtools' arguments.
OUTPUT_FORMAT = r"(--output-fmt|-O)+\s+(\S+)"

dragmap_hashtable_cmd = shell.create(
    name="dragmap_hashtable",
    image=DRAGMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"hashmap": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        ARGS={{inputs.args}}
        dragen-os --build-hash-table true --ht-reference "${{FA[0]}}" --output-directory {{outputs.hashmap}} \
            $ARGS --ht-num-threads $TASK_CPUS || exit 1
    """,
)

dragmap_align_cmd = shell.create(
    name="dragmap_align",
    image=DRAGMAP_SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "hashmap": Dir,
        "fasta": File | None,
        "prefix": str,
        "sort_bam": bool,
        "args": str,
        "args2": str,
    },
    defaults={"sort_bam": True, "args": "", "args2": ""},
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
        if [ ${{#R2[@]}} -eq 0 ]; then READS=(-1 "${{R1[0]}}"); else READS=(-1 "${{R1[0]}}" -2 "${{R2[0]}}"); fi
        SAMTOOLS=view
        if [ {{inputs.sort_bam}} = true ]; then SAMTOOLS=sort; fi

        dragen-os -r {{inputs.hashmap}} $ARGS --num-threads $TASK_CPUS "${{READS[@]}}" \
            2> "$OUT/$P.dragmap.log" \
            | samtools $SAMTOOLS $ARGS2 --threads $TASK_CPUS "${{REF[@]}}" -o "$OUT/$P.$EXT" -
        STATUS=$?
        cat "$OUT/$P.dragmap.log" >&2
        if [ $STATUS -ne 0 ]; then exit 1; fi
    """,
)

dragmap_hashtable_env = flyte.TaskEnvironment.from_task("dragmap_hashtable", dragmap_hashtable_cmd.as_task())
dragmap_align_env = flyte.TaskEnvironment.from_task("dragmap_align", dragmap_align_cmd.as_task())
env = flyte.TaskEnvironment(name="dragmap", depends_on=[dragmap_hashtable_env, dragmap_align_env])


async def dragmap_hashtable(fasta: File, args: str = "") -> Dir:
    """DRAGMAP's hash table of ``fasta``, as a Dir."""
    return await dragmap_hashtable_cmd(fasta=fasta, args=args)


@dataclass
class DragmapAlignResult:
    alignment: File  # <prefix>.bam, .sam or .cram
    alignment_index: File | None  # <prefix>.bam.csi or <prefix>.cram.crai, with --write-index in args2
    log: File  # <prefix>.dragmap.log
    results: Dir  # everything the task wrote


async def dragmap_align(
    reads_1: File,
    reads_2: File | None,
    hashmap: Dir,
    prefix: str,
    fasta: File | None = None,
    sort_bam: bool = True,
    args: str = "",
    args2: str = "",
) -> DragmapAlignResult:
    """Align one sample's reads with DRAGMAP; ``args`` go to DRAGMAP, ``args2`` to samtools.

    ``reads_2`` is None for single-end reads. ``fasta`` is the reference, needed
    for CRAM output (``--output-fmt cram`` in ``args2``).
    """
    found = re.search(OUTPUT_FORMAT, args2)
    extension = found.group(2).lower() if found else "bam"
    if extension == "cram" and fasta is None:
        raise ValueError("Fasta reference is required for CRAM output")
    results = await dragmap_align_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        hashmap=hashmap,
        fasta=fasta,
        prefix=prefix,
        sort_bam=sort_bam,
        args=args,
        args2=args2,
    )

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"dragmap align wrote no {name}")
        return file

    alignment_index = None
    for suffix in ("csi", "crai"):
        alignment_index = alignment_index or await results.get_file(f"{prefix}.{extension}.{suffix}")
    return DragmapAlignResult(
        alignment=await pick(f"{prefix}.{extension}"),
        alignment_index=alignment_index,
        log=await pick(f"{prefix}.dragmap.log"),
        results=results,
    )
