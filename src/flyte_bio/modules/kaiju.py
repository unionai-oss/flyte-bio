"""kaiju — protein-level taxonomic classification databases.

Exposes :func:`kaiju_mkfmi`, which builds a Kaiju index from a protein FASTA
(``kaiju-mkbwt``, then ``kaiju-mkfmi``). The database is a Dir holding, as
upstream's outputs, ``<prefix>.fmi`` with the taxonomy it goes with
(``nodes.dmp`` and ``names.dmp``, copied in), plus, with
``keep_intermediate``, the Burrows-Wheeler transform ``<prefix>.bwt`` and
suffix array ``<prefix>.sa``.

``-n`` is the task's CPUs; options such as the alphabet (``-a``) go in
``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KAIJU_IMAGE = "quay.io/biocontainers/kaiju:1.10.0--h43eeafb_0"
# Upstream labels it process_high; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

kaiju_mkfmi_cmd = shell.create(
    name="kaiju_mkfmi",
    image=KAIJU_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "nodes_dmp": File,
        "names_dmp": File,
        "keep_intermediate": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"keep_intermediate": False, "args": ""},
    outputs={"db": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}}); NODES=({{inputs.nodes_dmp}}); NAMES=({{inputs.names_dmp}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.db}}
        W=$(mktemp -d); cd "$W"
        kaiju-mkbwt \
            $ARGS \
            -n $TASK_CPUS \
            -o "$P" \
            "${{FASTA[0]}}" \
            || exit 1
        kaiju-mkfmi "$P" || exit 1
        cp "${{NODES[0]}}" "$OUT/nodes.dmp"
        cp "${{NAMES[0]}}" "$OUT/names.dmp"
        mv "$P.fmi" "$OUT"/
        if [ {{inputs.keep_intermediate}} = true ]; then mv "$P.bwt" "$P.sa" "$OUT"/; fi
        true
    """,
)

env = flyte.TaskEnvironment.from_task("kaiju", kaiju_mkfmi_cmd.as_task())


async def kaiju_mkfmi(
    fasta: File,
    nodes_dmp: File,
    names_dmp: File,
    keep_intermediate: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """Build a Kaiju index of protein ``fasta``: a Dir of ``<prefix>.fmi``, ``nodes.dmp`` and ``names.dmp``."""
    return await kaiju_mkfmi_cmd(
        fasta=fasta,
        nodes_dmp=nodes_dmp,
        names_dmp=names_dmp,
        keep_intermediate=keep_intermediate,
        prefix=prefix,
        args=args,
    )
