"""mafft — multiple sequence alignment with MAFFT.

Exposes :func:`mafft_align` and :func:`mafft_guidetree`.

:func:`mafft_align` aligns the sequences in ``fasta`` and writes the
alignment to ``<prefix>.fas``, or ``<prefix>.fas.gz`` (compressed with
pigz) with ``compress``. Any of the optional inputs adds sequences to the
alignment with the MAFFT option of the same name: ``add`` (``--add``),
``add_fragments`` (``--addfragments``), ``add_full`` (``--addfull``),
``add_profile`` (``--addprofile``) and ``add_long`` (``--addlong``).

As upstream, every input may be gzipped (each is read through ``unpigz -cdf``,
which passes uncompressed files through unchanged), and an input named
``<prefix>.fas`` is an error. MAFFT's own options (``--auto``, ``--maxiterate``)
go in ``args``.

:func:`mafft_guidetree` writes the guide tree MAFFT builds for ``fasta``
(``--retree 0 --treeout``) to ``<prefix>.dnd``. As upstream, the numeric
``<n>_`` prefixes MAFFT adds to the leaf names are removed, so other tools can
read the tree.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker image: mafft 7.520 and pigz.
MAFFT_IMAGE = (
    "quay.io/biocontainers/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:"
    "425707898cf4f85051b77848be253b88f1d2298a-0"
)
# Upstream's guide-tree image is a newer MAFFT (7.525) without pigz.
GUIDETREE_IMAGE = "quay.io/biocontainers/mafft:7.525--h031d066_1"
# Upstream labels this process_high; large alignments may need more.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

mafft_align_cmd = shell.create(
    name="mafft_align",
    image=MAFFT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "add": File | None,
        "add_fragments": File | None,
        "add_full": File | None,
        "add_profile": File | None,
        "add_long": File | None,
        "prefix": str,
        "compress": bool,
        "args": str,
    },
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        FA=({{inputs.fasta}})
        shopt -s nullglob
        ADD=({{inputs.add}}); ADD_FRAGMENTS=({{inputs.add_fragments}}); ADD_FULL=({{inputs.add_full}})
        ADD_PROFILE=({{inputs.add_profile}}); ADD_LONG=({{inputs.add_long}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        if [ "$(basename "${{FA[0]}}")" = "$P.fas" ]; then
            echo "Input and output names are the same, set prefix in module configuration to disambiguate!" >&2
            exit 1
        fi
        # Every input may be gzipped. MAFFT reads the alignment input from a
        # process substitution, as upstream; the added sequences are
        # decompressed to files first, since a process substitution in an
        # array assignment is closed before MAFFT runs. MAFFT can't write to a
        # process substitution, so the output is piped (with pipefail, unlike
        # upstream, so a MAFFT failure fails the task).
        EXTRA=()
        add_sequences() {{
            if [ $# -gt 1 ]; then
                unpigz -cdf "$2" > "$1.fa"
                EXTRA+=(--$1 "$1.fa")
            fi
        }}
        add_sequences add "${{ADD[@]}}"
        add_sequences addfragments "${{ADD_FRAGMENTS[@]}}"
        add_sequences addfull "${{ADD_FULL[@]}}"
        add_sequences addprofile "${{ADD_PROFILE[@]}}"
        add_sequences addlong "${{ADD_LONG[@]}}"

        if [ {{inputs.compress}} = true ]; then
            mafft --thread $TASK_CPUS "${{EXTRA[@]}}" $ARGS <(unpigz -cdf "${{FA[0]}}") \
                | pigz -cp $TASK_CPUS > "$OUT/$P.fas.gz"
        else
            mafft --thread $TASK_CPUS "${{EXTRA[@]}}" $ARGS <(unpigz -cdf "${{FA[0]}}") > "$OUT/$P.fas"
        fi
    """,
)

mafft_guidetree_cmd = shell.create(
    name="mafft_guidetree",
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

        # MAFFT writes the tree next to its input (<input>.tree), so the input
        # is linked into the working directory.
        NAME=$(basename "${{FA[0]}}")
        ln -s "${{FA[0]}}" "$NAME"
        mafft --thread $TASK_CPUS --retree 0 $ARGS --treeout "$NAME" > log.txt || exit 1
        mv ./*.tree "$P.dnd.tmp"

        # remove all prefixes added by mafft which make the output incompatible with other tools
        awk '{{gsub(/^[0-9]+_/, ""); print}}' "$P.dnd.tmp" > "$OUT/$P.dnd"
    """,
)

# One image per environment: the guide tree's MAFFT differs from the aligner's.
mafft_align_env = flyte.TaskEnvironment.from_task("mafft_align", mafft_align_cmd.as_task())
mafft_guidetree_env = flyte.TaskEnvironment.from_task("mafft_guidetree", mafft_guidetree_cmd.as_task())
env = flyte.TaskEnvironment(name="mafft", depends_on=[mafft_align_env, mafft_guidetree_env])


async def mafft_align(
    fasta: File,
    add: File | None = None,
    add_fragments: File | None = None,
    add_full: File | None = None,
    add_profile: File | None = None,
    add_long: File | None = None,
    prefix: str = "sample",
    compress: bool = False,
    args: str = "",
) -> File:
    """Align ``fasta`` (adding any ``add*`` sequences): ``<prefix>.fas``, or ``<prefix>.fas.gz`` with ``compress``."""
    results = await mafft_align_cmd(
        fasta=fasta,
        add=add,
        add_fragments=add_fragments,
        add_full=add_full,
        add_profile=add_profile,
        add_long=add_long,
        prefix=prefix,
        compress=compress,
        args=args,
    )
    name = f"{prefix}.fas.gz" if compress else f"{prefix}.fas"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"mafft wrote no {name}")
    return found


async def mafft_guidetree(fasta: File, prefix: str = "sample", args: str = "") -> File:
    """The guide tree MAFFT builds for ``fasta``, its leaf names unprefixed: ``<prefix>.dnd``."""
    results = await mafft_guidetree_cmd(fasta=fasta, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.dnd")
    if found is None:
        raise FileNotFoundError(f"mafft wrote no {prefix}.dnd")
    return found
