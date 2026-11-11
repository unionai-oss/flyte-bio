"""agat — GFF/GTF conversions with AGAT (Another Gtf/Gff Analysis Toolkit).

Exposes:

- :func:`agat_convertgff2bed`, which converts a GFF/GTF to BED12,
  ``<prefix>.bed`` (``agat_convert_sp_gff2bed.pl``).
- :func:`agat_convertspgff2tsv`, which converts a GFF/GTF to a table of one
  feature per row, ``<prefix>.tsv`` (``agat_convert_sp_gff2tsv.pl``).
- :func:`agat_spkeeplongestisoform`, which keeps each gene's longest isoform,
  ``<prefix>.longest.gff`` (``agat_sp_keep_longest_isoform.pl``), optionally
  with an AGAT config YAML. ``--cpu`` follows the task's CPUs.
- :func:`agat_convertspgxf2gxf`, which checks, fixes and standardises a
  GFF/GTF into GFF3, ``<prefix>.agat.gff`` (``agat_convert_sp_gxf2gxf.pl``).
- :func:`agat_spstatistics`, which summarises an annotation's features,
  ``<prefix>.stats.txt`` (``agat_sp_statistics.pl``).

All share one image. AGAT writes its log (``<input>.agat.log``) to the working
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

# Upstream's process_medium; AGAT parallelises over its --cpu.
KEEPLONGEST_RESOURCES = flyte.Resources(cpu=4, memory="12Gi")

agat_spkeeplongestisoform_cmd = shell.create(
    name="agat_spkeeplongestisoform",
    image=AGAT_IMAGE,
    resources=KEEPLONGEST_RESOURCES,
    inputs={"gxf": File, "config": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GXF=({{inputs.gxf}})
        shopt -s nullglob; CONFIG=({{inputs.config}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{GXF[0]}}" "$(basename "${{GXF[0]}}")"
        CONFIG_ARGS=()
        if [ ${{#CONFIG[@]}} -gt 0 ]; then
            ln -s "${{CONFIG[0]}}" "$(basename "${{CONFIG[0]}}")"
            CONFIG_ARGS=(--config "$(basename "${{CONFIG[0]}}")")
        fi
        agat_sp_keep_longest_isoform.pl \
            --gff "$(basename "${{GXF[0]}}")" \
            --cpu $TASK_CPUS \
            "${{CONFIG_ARGS[@]}}" \
            --out "$P.longest.gff" \
            $ARGS || exit 1
        mv "$P.longest.gff" {{outputs.results}}/
    """,
)

agat_convertspgxf2gxf_cmd = shell.create(
    name="agat_convertspgxf2gxf",
    image=AGAT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gxf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GXF=({{inputs.gxf}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{GXF[0]}}" "$(basename "${{GXF[0]}}")"
        agat_convert_sp_gxf2gxf.pl \
            --gxf "$(basename "${{GXF[0]}}")" \
            --output "$P.agat.gff" \
            $ARGS || exit 1
        mv "$P.agat.gff" {{outputs.results}}/
    """,
)

# Upstream's process_low.
STATISTICS_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")

agat_spstatistics_cmd = shell.create(
    name="agat_spstatistics",
    image=AGAT_IMAGE,
    resources=STATISTICS_RESOURCES,
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
        agat_sp_statistics.pl \
            --gff "$(basename "${{GFF[0]}}")" \
            --output "$P.stats.txt" \
            $ARGS || exit 1
        mv "$P.stats.txt" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "agat",
    agat_convertgff2bed_cmd.as_task(),
    agat_convertspgff2tsv_cmd.as_task(),
    agat_spkeeplongestisoform_cmd.as_task(),
    agat_convertspgxf2gxf_cmd.as_task(),
    agat_spstatistics_cmd.as_task(),
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


async def agat_spkeeplongestisoform(gxf: File, prefix: str, config: File | None = None, args: str = "") -> File:
    """``gxf`` with only each gene's longest isoform: ``<prefix>.longest.gff``."""
    if gxf.name == f"{prefix}.longest.gff":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await agat_spkeeplongestisoform_cmd(gxf=gxf, config=config, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.longest.gff", "sp_keep_longest_isoform")


async def agat_convertspgxf2gxf(gxf: File, prefix: str, args: str = "") -> File:
    """``gxf`` checked, fixed and written as GFF3: ``<prefix>.agat.gff``."""
    if gxf.name == f"{prefix}.agat.gff":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await agat_convertspgxf2gxf_cmd(gxf=gxf, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.agat.gff", "convert_sp_gxf2gxf")


async def agat_spstatistics(gff: File, prefix: str, args: str = "") -> File:
    """Feature statistics of ``gff``: ``<prefix>.stats.txt``."""
    results = await agat_spstatistics_cmd(gff=gff, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.stats.txt", "sp_statistics")
