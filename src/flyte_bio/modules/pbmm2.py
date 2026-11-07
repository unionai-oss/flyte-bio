"""pbmm2 — PacBio's minimap2 wrapper: alignment of PacBio reads (BAM, FASTQ, FASTA) to a reference.

Exposes :func:`pbmm2_align`, which runs ``pbmm2 align`` on one sample's reads
against a reference FASTA (or ``.mmi`` index) and writes ``<prefix>.bam``.
With ``--sort`` and ``--bam-index BAI`` (or ``CSI``) in ``args``, pbmm2 also
indexes it, to ``<prefix>.bam.bai`` (or ``.csi``).

As upstream:

- A ``.fna`` or ``.fna.gz`` reference is renamed to ``.fa`` / ``.fa.gz``,
  which pbmm2 recognises as FASTA.
- pbmm2 gets the task's CPUs as ``--num-threads``.
- Inputs are linked into a working directory under their own names, so the
  BAM header's ``@PG`` command line names them as upstream's do.

``args`` go to ``pbmm2 align``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PBMM2_IMAGE = "community.wave.seqera.io/library/pbmm2:26.2.0--37598eea709c00f6"
# Upstream labels it process_high; its tests ran on 4 CPUs, which the output
# BAM's @PG command line records (--num-threads). A human index needs ~16 GB.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

pbmm2_align_cmd = shell.create(
    name="pbmm2_align",
    image=PBMM2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "reference": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}}); REF=({{inputs.reference}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{READS[0]}}"); ln -s "${{READS[0]}}" "$NAME"
        REF_NAME=$(basename "${{REF[0]}}")
        case "$REF_NAME" in
            *.fna) REF_NAME="${{REF_NAME%.fna}}.fa" ;;
            *.fna.gz) REF_NAME="${{REF_NAME%.fna.gz}}.fa.gz" ;;
        esac
        ln -s "${{REF[0]}}" "$REF_NAME"
        pbmm2 align $ARGS "$REF_NAME" "$NAME" "$P.bam" --num-threads $TASK_CPUS || exit 1
        for f in "$P.bam" "$P.bam.bai" "$P.bam.csi"; do
            if [ -f "$f" ]; then mv "$f" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("pbmm2", pbmm2_align_cmd.as_task())


@dataclass
class Pbmm2AlignResult:
    bam: File  # <prefix>.bam
    index: File | None  # <prefix>.bam.bai or .csi, with --sort and --bam-index
    results: Dir  # both


async def pbmm2_align(reads: File, reference: File, prefix: str = "sample", args: str = "") -> Pbmm2AlignResult:
    """Align ``reads`` (PacBio BAM, FASTQ or FASTA) to ``reference`` (FASTA or ``.mmi``): ``<prefix>.bam``."""
    results = await pbmm2_align_cmd(reads=reads, reference=reference, prefix=prefix, args=args)
    bam = await results.get_file(f"{prefix}.bam")
    if bam is None:
        raise FileNotFoundError(f"pbmm2 align wrote no {prefix}.bam")
    index = await results.get_file(f"{prefix}.bam.bai") or await results.get_file(f"{prefix}.bam.csi")
    return Pbmm2AlignResult(bam=bam, index=index, results=results)
