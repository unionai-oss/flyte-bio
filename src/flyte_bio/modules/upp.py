"""upp — multiple sequence alignment with UPP (Ultra-large alignments using Phylogeny-aware Profiles).

Exposes :func:`upp_align`, which aligns ``fasta_unaligned``, and/or extends a
backbone alignment ``fasta_aligned`` (``--alignment``) with it, and writes
``<prefix>.aln``, or ``<prefix>.aln.gz`` (gzipped) with ``compress``. With
``tree`` (Newick), UPP uses it as the backbone tree (``-t``). UPP runs with
``-x`` set to the task's CPUs; its other options (``-m amino``) go in
``args``. Its temporary and working files go to a scratch directory.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

UPP_IMAGE = "quay.io/biocontainers/sepp:4.5.6--py312h87e0c26_4"
# Upstream labels this process_medium; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

upp_align_cmd = shell.create(
    name="upp_align",
    image=UPP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta_unaligned": File | None,
        "fasta_aligned": File | None,
        "tree": File | None,
        "prefix": str,
        "compress": bool,
        "args": str,
    },
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        UNALIGNED=({{inputs.fasta_unaligned}}); ALIGNED=({{inputs.fasta_aligned}}); TREE=({{inputs.tree}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        TREE_ARGS=(); SEQ_CMD=(); ALIGN_CMD=()
        if [ ${{#TREE[@]}} -gt 0 ]; then TREE_ARGS=(-t "${{TREE[0]}}"); fi
        if [ ${{#UNALIGNED[@]}} -gt 0 ]; then SEQ_CMD=(--sequence_file "${{UNALIGNED[0]}}"); fi
        if [ ${{#ALIGNED[@]}} -gt 0 ]; then ALIGN_CMD=(--alignment "${{ALIGNED[0]}}"); fi

        run_upp.py \
            $ARGS \
            "${{TREE_ARGS[@]}}" \
            "${{SEQ_CMD[@]}}" \
            "${{ALIGN_CMD[@]}}" \
            -x $TASK_CPUS \
            -p ./upp_tmp/ \
            -d ./upp_output/ \
            -o "$P" || exit 1

        mv upp_output/"$P"_alignment.fasta "$P.aln"
        if [ {{inputs.compress}} = true ]; then
            gzip "$P.aln"
        fi
        mv "$P".aln* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("upp", upp_align_cmd.as_task())


async def upp_align(
    fasta_unaligned: File | None = None,
    fasta_aligned: File | None = None,
    tree: File | None = None,
    prefix: str = "sample",
    compress: bool = False,
    args: str = "",
) -> File:
    """Align ``fasta_unaligned`` (onto the backbone ``fasta_aligned``/``tree``, if given): ``<prefix>.aln(.gz)``."""
    if fasta_unaligned is None and fasta_aligned is None:
        raise ValueError("upp_align needs fasta_unaligned, fasta_aligned or both")
    results = await upp_align_cmd(
        fasta_unaligned=fasta_unaligned,
        fasta_aligned=fasta_aligned,
        tree=tree,
        prefix=prefix,
        compress=compress,
        args=args,
    )
    name = f"{prefix}.aln.gz" if compress else f"{prefix}.aln"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"upp wrote no {name}")
    return found
