"""clustalo — multiple sequence alignment with Clustal Omega.

Exposes:

- :func:`clustalo_align` — align the sequences in ``fasta``, writing
  ``<prefix>.aln``, or ``<prefix>.aln.gz`` (compressed with pigz) with
  ``compress``. The optional inputs map to Clustal Omega's options of the
  same name: ``tree`` (``--guidetree-in``), ``hmm_in`` (``--hmm-in``),
  ``hmm_batch`` (``--hmm-batch``), ``profile1`` (``--profile1``) and
  ``profile2`` (``--profile2``). The output format (``--outfmt``) and other
  options go in ``args``.
- :func:`clustalo_guidetree` — the guide tree Clustal Omega builds for
  ``fasta``, ``<prefix>.dnd`` (Newick), for :func:`clustalo_align`'s ``tree``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker images: Clustal Omega 1.2.4 with pigz for alignments, and
# Clustal Omega alone for guide trees.
CLUSTALO_IMAGE = (
    "quay.io/biocontainers/mulled-v2-4cefc38542f86c17596c29b35a059de10387c6a7:"
    "adbe4fbad680f9beb083956d79128039a727e7b3-0"
)
GUIDETREE_IMAGE = "quay.io/biocontainers/clustalo:1.2.4--h87f3376_5"
# Upstream labels both process_medium; large alignments may need more.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

clustalo_align_cmd = shell.create(
    name="clustalo_align",
    image=CLUSTALO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "tree": File | None,
        "hmm_in": File | None,
        "hmm_batch": File | None,
        "profile1": File | None,
        "profile2": File | None,
        "prefix": str,
        "compress": bool,
        "args": str,
    },
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob
        TREE=({{inputs.tree}}); HMM_IN=({{inputs.hmm_in}}); HMM_BATCH=({{inputs.hmm_batch}})
        PROFILE1=({{inputs.profile1}}); PROFILE2=({{inputs.profile2}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        EXTRA=()
        if [ ${{#TREE[@]}} -gt 0 ]; then EXTRA+=("--guidetree-in=${{TREE[0]}}"); fi
        if [ ${{#HMM_IN[@]}} -gt 0 ]; then EXTRA+=("--hmm-in=${{HMM_IN[0]}}"); fi
        if [ ${{#HMM_BATCH[@]}} -gt 0 ]; then EXTRA+=("--hmm-batch=${{HMM_BATCH[0]}}"); fi
        if [ ${{#PROFILE1[@]}} -gt 0 ]; then EXTRA+=("--profile1=${{PROFILE1[0]}}"); fi
        if [ ${{#PROFILE2[@]}} -gt 0 ]; then EXTRA+=("--profile2=${{PROFILE2[0]}}"); fi

        if [ {{inputs.compress}} = true ]; then
            # As upstream: the shell doesn't wait for a process substitution, so
            # wait for pigz; --force lets Clustal Omega open the existing
            # /dev/fd path.
            clustalo -i "${{FA[0]}}" "${{EXTRA[@]}}" --threads=$TASK_CPUS $ARGS \
                --force -o >(pigz -cp $TASK_CPUS > "$OUT/$P.aln.gz") && wait $!
        else
            clustalo -i "${{FA[0]}}" "${{EXTRA[@]}}" --threads=$TASK_CPUS $ARGS -o "$OUT/$P.aln"
        fi
    """,
)

clustalo_guidetree_cmd = shell.create(
    name="clustalo_guidetree",
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

        clustalo -i "${{FA[0]}}" --guidetree-out "$OUT/$P.dnd" --threads=$TASK_CPUS $ARGS
    """,
)

# One image per environment.
clustalo_align_env = flyte.TaskEnvironment.from_task("clustalo_align", clustalo_align_cmd.as_task())
clustalo_guidetree_env = flyte.TaskEnvironment.from_task("clustalo_guidetree", clustalo_guidetree_cmd.as_task())
env = flyte.TaskEnvironment(name="clustalo", depends_on=[clustalo_align_env, clustalo_guidetree_env])


async def clustalo_align(
    fasta: File,
    tree: File | None = None,
    hmm_in: File | None = None,
    hmm_batch: File | None = None,
    profile1: File | None = None,
    profile2: File | None = None,
    prefix: str = "sample",
    compress: bool = False,
    args: str = "",
) -> File:
    """Align ``fasta``: ``<prefix>.aln``, or ``<prefix>.aln.gz`` with ``compress``."""
    results = await clustalo_align_cmd(
        fasta=fasta,
        tree=tree,
        hmm_in=hmm_in,
        hmm_batch=hmm_batch,
        profile1=profile1,
        profile2=profile2,
        prefix=prefix,
        compress=compress,
        args=args,
    )
    name = f"{prefix}.aln.gz" if compress else f"{prefix}.aln"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"clustalo wrote no {name}")
    return found


async def clustalo_guidetree(fasta: File, prefix: str = "sample", args: str = "") -> File:
    """The guide tree for ``fasta``: ``<prefix>.dnd``."""
    results = await clustalo_guidetree_cmd(fasta=fasta, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.dnd")
    if found is None:
        raise FileNotFoundError(f"clustalo wrote no {prefix}.dnd")
    return found
