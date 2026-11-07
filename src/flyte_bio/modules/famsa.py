"""famsa — fast multiple sequence alignment with FAMSA.

Exposes:

- :func:`famsa_align` — align the sequences in ``fasta``, writing
  ``<prefix>.aln``, or ``<prefix>.aln.gz`` (FAMSA's own ``-gz``) with
  ``compress``. With ``tree`` (a Newick guide tree, e.g. from
  :func:`famsa_guidetree`), FAMSA imports it (``-gt import``) instead of
  building its own.
- :func:`famsa_guidetree` — the guide tree FAMSA builds for ``fasta``,
  ``<prefix>.dnd`` (Newick, ``-gt_export``).

FAMSA's options (``-gt upgma``, ``-medoidtree``) go in ``args``. Its output
depends on the thread count (``-t``, the task's CPUs): the defaults' 4 CPUs
match upstream's tests.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FAMSA_IMAGE = "quay.io/biocontainers/famsa:2.4.1--h9ee0642_0"
# Upstream labels both process_medium. The alignment depends on the thread
# count; 4 CPUs is what upstream's tests run with.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

famsa_align_cmd = shell.create(
    name="famsa_align",
    image=FAMSA_IMAGE,
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

        EXTRA=()
        if [ ${{#TREE[@]}} -gt 0 ]; then EXTRA+=(-gt import "${{TREE[0]}}"); fi
        NAME="$P.aln"
        if [ {{inputs.compress}} = true ]; then EXTRA+=(-gz); NAME="$P.aln.gz"; fi
        famsa "${{EXTRA[@]}}" $ARGS -t $TASK_CPUS "${{FA[0]}}" "$OUT/$NAME" || exit 1
    """,
)

famsa_guidetree_cmd = shell.create(
    name="famsa_guidetree",
    image=FAMSA_IMAGE,
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

        famsa -gt_export $ARGS -t $TASK_CPUS "${{FA[0]}}" "$OUT/$P.dnd" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("famsa", famsa_align_cmd.as_task(), famsa_guidetree_cmd.as_task())


async def famsa_align(
    fasta: File, tree: File | None = None, prefix: str = "sample", compress: bool = False, args: str = ""
) -> File:
    """Align ``fasta`` (with ``tree`` as the guide tree, if given): ``<prefix>.aln`` or ``<prefix>.aln.gz``."""
    results = await famsa_align_cmd(fasta=fasta, tree=tree, prefix=prefix, compress=compress, args=args)
    name = f"{prefix}.aln.gz" if compress else f"{prefix}.aln"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"famsa wrote no {name}")
    return found


async def famsa_guidetree(fasta: File, prefix: str = "sample", args: str = "") -> File:
    """The guide tree FAMSA builds for ``fasta``: ``<prefix>.dnd``."""
    results = await famsa_guidetree_cmd(fasta=fasta, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.dnd")
    if found is None:
        raise FileNotFoundError(f"famsa wrote no {prefix}.dnd")
    return found
