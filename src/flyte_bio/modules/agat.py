"""agat — GFF/GTF conversions with AGAT (Another Gtf/Gff Analysis Toolkit).

Exposes:

- :func:`agat_convertgff2bed`, which converts a GFF/GTF to BED12,
  ``<prefix>.bed`` (``agat_convert_sp_gff2bed.pl``).
- :func:`agat_convertspgff2tsv`, which converts a GFF/GTF to a table of one
  feature per row, ``<prefix>.tsv`` (``agat_convert_sp_gff2tsv.pl``).

Both share one image. AGAT writes its log (``<input>.agat.log``) to the working
directory, which is not kept, as upstream doesn't emit it.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

AGAT_IMAGE = "community.wave.seqera.io/library/agat:1.6.1--1f8d94d325805bb4"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

agat_convertgff2bed_cmd = shell.create(
    name="agat_convertgff2bed",
    image=AGAT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GFF=({{inputs.gff}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{GFF[0]}}" "$(basename "${{GFF[0]}}")"
        agat_convert_sp_gff2bed.pl \
            --gff "$(basename "${{GFF[0]}}")" \
            -o "$P.bed" \
            $ARGS || exit 1
        mv "$P.bed" {{outputs.results}}/
    """,
)

agat_convertspgff2tsv_cmd = shell.create(
    name="agat_convertspgff2tsv",
    image=AGAT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GFF=({{inputs.gff}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{GFF[0]}}" "$(basename "${{GFF[0]}}")"
        agat_convert_sp_gff2tsv.pl \
            --gff "$(basename "${{GFF[0]}}")" \
            --output "$P.tsv" \
            $ARGS || exit 1
        mv "$P.tsv" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "agat",
    agat_convertgff2bed_cmd.as_task(),
    agat_convertspgff2tsv_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"agat {tool} wrote no {name}")
    return found


async def agat_convertgff2bed(gff: File, prefix: str, args: str = "") -> File:
    """``gff`` as BED12: ``<prefix>.bed``. Feature types to keep (``-f gene,mRNA,…``) go in ``args``."""
    if gff.name == f"{prefix}.bed":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await agat_convertgff2bed_cmd(gff=gff, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bed", "convert_sp_gff2bed")


async def agat_convertspgff2tsv(gff: File, prefix: str, args: str = "") -> File:
    """``gff`` as a table of one feature per row: ``<prefix>.tsv``."""
    if gff.name == f"{prefix}.tsv":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await agat_convertspgff2tsv_cmd(gff=gff, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.tsv", "convert_sp_gff2tsv")
