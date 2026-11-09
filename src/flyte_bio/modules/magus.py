"""magus — multiple sequence alignment with MAGUS.

Exposes:

- :func:`magus_align` — align the sequences in ``fasta``, writing
  ``<prefix>.aln``, or ``<prefix>.aln.gz`` (gzipped with pigz) with
  ``compress``. With ``tree`` (Newick, e.g. from :func:`magus_guidetree`),
  MAGUS uses it to decompose the sequences (``-t``).
- :func:`magus_guidetree` — the guide tree MAGUS builds for ``fasta``
  (``--onlyguidetree TRUE``): ``<prefix>.tree``.

Both run MAGUS with ``-np`` set to the task's CPUs; its other options go in
``args``. MAGUS's working files go to a scratch directory. Upstream writes the
compressed alignment through a process substitution (with ``--overwrite``);
here MAGUS writes the file and pigz compresses it after, which is the same
alignment. MAGUS's output differs between machines (upstream's tests check
only that every sequence is there).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker images: MAGUS 0.2.0 with pigz for alignments, and MAGUS
# alone for guide trees.
MAGUS_IMAGE = (
    "quay.io/biocontainers/mulled-v2-ae4ea1182e75371808710b6c081bef8b228c4815:"
    "10b41722a6b9471a0945fe6baeb9aff444d8eb1d-0"
)
GUIDETREE_IMAGE = "quay.io/biocontainers/magus-msa:0.2.0--pyhdfd78af_0"
# Upstream: align process_high, guidetree process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

magus_align_cmd = shell.create(
    name="magus_align",
    image=MAGUS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "tree": File | None, "prefix": str, "compress": bool, "args": str},
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; TREE=({{inputs.tree}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        LOADTREE=()
        if [ ${{#TREE[@]}} -gt 0 ]; then LOADTREE=(-t "${{TREE[0]}}"); fi
        magus \
            -np $TASK_CPUS \
            -i "${{FA[0]}}" \
            -d ./ \
            -o "$P.aln" \
            "${{LOADTREE[@]}}" \
            $ARGS || exit 1
        if [ {{inputs.compress}} = true ]; then
            pigz -p $TASK_CPUS "$P.aln"
        fi
        mv "$P".aln* "$OUT"/
    """,
)

magus_guidetree_cmd = shell.create(
    name="magus_guidetree",
    image=GUIDETREE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        magus \
            -np $TASK_CPUS \
            -i "${{FA[0]}}" \
            -o "$OUT/$P.tree" \
            --onlyguidetree TRUE \
            $ARGS || exit 1
    """,
)

# One image per environment: the aligner's adds pigz.
magus_align_env = flyte.TaskEnvironment.from_task("magus_align", magus_align_cmd.as_task())
magus_guidetree_env = flyte.TaskEnvironment.from_task("magus_guidetree", magus_guidetree_cmd.as_task())
env = flyte.TaskEnvironment(name="magus", depends_on=[magus_align_env, magus_guidetree_env])


async def magus_align(
    fasta: File, tree: File | None = None, prefix: str = "sample", compress: bool = False, args: str = ""
) -> File:
    """Align ``fasta`` (decomposed along ``tree``, if given): ``<prefix>.aln`` or ``<prefix>.aln.gz``."""
    results = await magus_align_cmd(fasta=fasta, tree=tree, prefix=prefix, compress=compress, args=args)
    name = f"{prefix}.aln.gz" if compress else f"{prefix}.aln"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"magus wrote no {name}")
    return found


async def magus_guidetree(fasta: File, prefix: str = "sample", args: str = "") -> File:
    """The guide tree MAGUS builds for ``fasta``: ``<prefix>.tree``."""
    results = await magus_guidetree_cmd(fasta=fasta, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.tree")
    if found is None:
        raise FileNotFoundError(f"magus wrote no {prefix}.tree")
    return found
