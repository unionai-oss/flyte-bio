"""muscle5 — multiple sequence alignment with MUSCLE 5's Super5 algorithm.

Exposes :func:`muscle5_super5`, which aligns the sequences in ``fasta`` with
``muscle -super5`` and writes ``<prefix>.aln``, or ``<prefix>.aln.gz``
(gzipped with pigz) with ``compress``. MUSCLE runs with ``-threads`` set to the
task's CPUs; its other options go in ``args``.

With ``-perm all`` in ``args``, MUSCLE writes one alignment per guide-tree
permutation: as upstream, the output name becomes ``<prefix>@.aln``, which
MUSCLE expands to ``<prefix><perm>.<n>.aln`` (``<prefix>abc.0.aln``,
``<prefix>none.0.aln``, …), each gzipped on its own with ``compress``.

Upstream writes a single compressed alignment through a process substitution;
here MUSCLE writes the file and pigz compresses it after, which is the same
alignment.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker image: MUSCLE 5.1 and pigz.
MUSCLE5_IMAGE = (
    "quay.io/biocontainers/mulled-v2-8eb01a3c2755c935d070dd03ff2dee698eeb4466:"
    "ceb6e65e00346ed20d0d8078dddf9858a7af0fe2-0"
)
# Upstream labels this process_medium; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

muscle5_super5_cmd = shell.create(
    name="muscle5_super5",
    image=MUSCLE5_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "compress": bool, "args": str},
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        if [[ " $ARGS " == *" -perm all "* ]]; then P="$P@"; fi
        # -output stays last, as upstream
        muscle \
            -super5 "${{FA[0]}}" \
            $ARGS \
            -threads $TASK_CPUS \
            -output "$P.aln" || exit 1

        # output may be multiple files if -perm all is set
        # compress these individually if set to compress output
        if [ {{inputs.compress}} = true ]; then
            pigz -p $TASK_CPUS ./*.aln
        fi
        mv ./*.aln* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("muscle5", muscle5_super5_cmd.as_task())


async def muscle5_super5(fasta: File, prefix: str = "sample", compress: bool = False, args: str = "") -> list[File]:
    """Align ``fasta`` with Super5: ``[<prefix>.aln(.gz)]``, or one alignment per permutation with ``-perm all``."""
    results = await muscle5_super5_cmd(fasta=fasta, prefix=prefix, compress=compress, args=args)
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)
    if not files:
        raise FileNotFoundError("muscle wrote no alignment")
    return files
