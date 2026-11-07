"""cmaple — maximum-likelihood phylogenies of closely related sequences with CMAPLE.

Exposes :func:`cmaple`, which infers a tree from the alignment ``aln``
(optionally gzipped, as ``.gz``; decompressed first, as upstream), starting
from the Newick tree ``newick`` (``-t``) if given. As upstream it runs the
``cmaple-aa`` build, which handles DNA as well as protein alignments.
Outputs keep upstream's names: ``<prefix>.treefile`` (Newick) and CMAPLE's
log, ``<prefix>.log`` (it records the command line and run times).
CMAPLE's options (``-st``, ``-m``) go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CMAPLE_IMAGE = "quay.io/biocontainers/cmaple:1.1.0--h503566f_1"
# Upstream labels this process_medium; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

cmaple_cmd = shell.create(
    name="cmaple",
    image=CMAPLE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"aln": File, "newick": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ALN=({{inputs.aln}})
        shopt -s nullglob; TREE=({{inputs.newick}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        NAME=$(basename "${{ALN[0]}}")
        if [ "${{NAME%.gz}}" != "$NAME" ]; then
            gzip -c -d "${{ALN[0]}}" > "${{NAME%.gz}}" || exit 1
            INPUT="${{NAME%.gz}}"
        else
            INPUT="${{ALN[0]}}"
        fi
        TREEOPT=()
        if [ ${{#TREE[@]}} -gt 0 ]; then TREEOPT=(-t "${{TREE[0]}}"); fi
        cmaple-aa \
            $ARGS \
            -nt $TASK_CPUS \
            --prefix "$P" \
            "${{TREEOPT[@]}}" \
            -aln "$INPUT" || exit 1
        mv "$P.treefile" "$P.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("cmaple", cmaple_cmd.as_task())


@dataclass
class CmapleResult:
    treefile: File  # <prefix>.treefile, the tree (Newick)
    log: File  # <prefix>.log


async def cmaple(aln: File, newick: File | None = None, prefix: str = "sample", args: str = "") -> CmapleResult:
    """The maximum-likelihood tree of ``aln`` (from ``newick``, if given), with CMAPLE's log."""
    results = await cmaple_cmd(aln=aln, newick=newick, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"cmaple wrote no {name}")
        return found

    return CmapleResult(treefile=await pick(f"{prefix}.treefile"), log=await pick(f"{prefix}.log"))
