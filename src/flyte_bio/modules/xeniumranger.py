"""xeniumranger — 10x Genomics Xenium Ranger: off-instrument reanalysis of Xenium In Situ output bundles.

Exposes:

- :func:`xeniumranger_relabel` — change the gene labels applied to decoded
  transcripts, from a new gene panel (``xeniumranger relabel``).
- :func:`xeniumranger_resegment` — resegment cells with Xenium Ranger's own
  segmentation, tuned through ``args`` (``xeniumranger resegment``).
- :func:`xeniumranger_importsegmentation` — import a segmentation made by
  another tool, either image-based (``nuclei`` and/or ``cells`` masks or
  polygons) or transcript-based (``transcript_assignment`` and
  ``viz_polygons``), and reassign transcripts to it
  (``xeniumranger import-segmentation``).

A Xenium bundle holds many files and is often large, so bundles go in and come
out as ``.tar.gz`` archives rather than Dirs. An input archive holds the
bundle's files either at its top level or in one folder; the bundle is the
folder holding ``experiment.xenium``. Each tool returns Xenium Ranger's
``outs`` folder, renamed ``<prefix>``, as ``<prefix>.tar.gz``: itself a
bundle, so tools can be chained.

As upstream, Xenium Ranger's cores and memory (``--localcores``,
``--localmem``, in GB) are the task's. For import-segmentation, as upstream:

- image-based and transcript-based inputs are mutually exclusive;
- ``expansion_distance`` is passed only with ``nuclei`` (Xenium Ranger rejects
  it otherwise);
- a ``coordinate_transform`` forces ``--units=microns``; otherwise ``units``
  is passed when given.

Telemetry is disabled (``TENX_DISABLE_TELEMETRY``). The container is
upstream's build of Xenium Ranger 4.0, distributed under 10x Genomics'
end-user license agreement.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

XENIUMRANGER_IMAGE = "quay.io/nf-core/xeniumranger:4.0"
NO_TELEMETRY = {"TENX_DISABLE_TELEMETRY": "1"}
# Upstream's process_high asks for 12 CPUs and 72 GB; its test bundles run in
# a few CPUs and GiB. Raise for full-size bundles with <cmd>.override(resources=...).
XENIUMRANGER_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")

# Extracts the bundle archive and sets BUNDLE to the folder holding experiment.xenium.
EXTRACT_BUNDLE = r"""
        mkdir bundle_archive
        tar -xzf "${BUNDLE_ARCHIVE[0]}" -C bundle_archive || exit 1
        BUNDLE=$(find "$W/bundle_archive" -maxdepth 2 -name experiment.xenium -printf '%h\n' | head -n 1)
        if [ -z "$BUNDLE" ]; then echo "no experiment.xenium in the bundle archive" >&2; exit 1; fi
"""

# Packs Xenium Ranger's outs (run id $RUN_ID) as <prefix>.tar.gz.
PACK_OUTS = r"""
        rm -rf "$P"
        mv "$RUN_ID/outs" "$P"
        tar -czf "$OUT/$P.tar.gz" "$P"
"""

xeniumranger_relabel_cmd = shell.create(
    name="xeniumranger_relabel",
    image=XENIUMRANGER_IMAGE,
    resources=XENIUMRANGER_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"bundle": File, "panel": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BUNDLE_ARCHIVE=({{inputs.bundle}})
        PANEL=({{inputs.panel}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        {EXTRACT_BUNDLE}
        RUN_ID=XENIUMRANGER_RELABEL
        xeniumranger relabel \
            --id="$RUN_ID" \
            --xenium-bundle="$BUNDLE" \
            --panel="${{PANEL[0]}}" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS || exit 1
        {PACK_OUTS}
    """,
)

xeniumranger_resegment_cmd = shell.create(
    name="xeniumranger_resegment",
    image=XENIUMRANGER_IMAGE,
    resources=XENIUMRANGER_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"bundle": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BUNDLE_ARCHIVE=({{inputs.bundle}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        {EXTRACT_BUNDLE}
        RUN_ID=XENIUMRANGER_RESEGMENT
        xeniumranger resegment \
            --id="$RUN_ID" \
            --xenium-bundle="$BUNDLE" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS || exit 1
        {PACK_OUTS}
    """,
)

xeniumranger_importsegmentation_cmd = shell.create(
    name="xeniumranger_importsegmentation",
    image=XENIUMRANGER_IMAGE,
    resources=XENIUMRANGER_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={
        "bundle": File,
        "transcript_assignment": File | None,
        "viz_polygons": File | None,
        "nuclei": File | None,
        "cells": File | None,
        "coordinate_transform": File | None,
        "units": str,
        "expansion_distance": str,
        "prefix": str,
        "args": str,
    },
    defaults={"units": "", "expansion_distance": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BUNDLE_ARCHIVE=({{inputs.bundle}})
        shopt -s nullglob
        ASSIGNMENT=({{inputs.transcript_assignment}}); POLYGONS=({{inputs.viz_polygons}})
        NUCLEI=({{inputs.nuclei}}); CELLS=({{inputs.cells}}); TRANSFORM=({{inputs.coordinate_transform}})
        shopt -u nullglob
        UNITS={{inputs.units}}
        EXPANSION={{inputs.expansion_distance}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        {EXTRACT_BUNDLE}

        EXTRA=()
        if [ ${{#NUCLEI[@]}} -gt 0 ] && [ -n "$EXPANSION" ]; then EXTRA+=(--expansion-distance="$EXPANSION"); fi
        if [ ${{#NUCLEI[@]}} -gt 0 ]; then EXTRA+=(--nuclei="${{NUCLEI[0]}}"); fi
        if [ ${{#CELLS[@]}} -gt 0 ]; then EXTRA+=(--cells="${{CELLS[0]}}"); fi
        if [ ${{#ASSIGNMENT[@]}} -gt 0 ]; then EXTRA+=(--transcript-assignment="${{ASSIGNMENT[0]}}"); fi
        if [ ${{#POLYGONS[@]}} -gt 0 ]; then EXTRA+=(--viz-polygons="${{POLYGONS[0]}}"); fi
        if [ ${{#TRANSFORM[@]}} -gt 0 ]; then
            EXTRA+=(--coordinate-transform="${{TRANSFORM[0]}}" --units=microns)
        elif [ -n "$UNITS" ]; then
            EXTRA+=(--units="$UNITS")
        fi

        RUN_ID=XENIUMRANGER_IMPORTSEGMENTATION
        xeniumranger import-segmentation \
            --id="$RUN_ID" \
            --xenium-bundle="$BUNDLE" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS \
            "${{EXTRA[@]}}" || exit 1
        {PACK_OUTS}
    """,
)

env = flyte.TaskEnvironment.from_task(
    "xeniumranger",
    xeniumranger_relabel_cmd.as_task(),
    xeniumranger_resegment_cmd.as_task(),
    xeniumranger_importsegmentation_cmd.as_task(),
)


async def xeniumranger_outs(results: Dir, prefix: str, tool: str) -> File:
    archive = await results.get_file(f"{prefix}.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"xeniumranger {tool} wrote no {prefix}.tar.gz")
    return archive


async def xeniumranger_relabel(bundle: File, panel: File, prefix: str = "xenium", args: str = "") -> File:
    """``bundle`` (a ``.tar.gz``) with its transcripts relabelled from gene panel ``panel``: ``<prefix>.tar.gz``."""
    results = await xeniumranger_relabel_cmd(bundle=bundle, panel=panel, prefix=prefix, args=args)
    return await xeniumranger_outs(results, prefix, "relabel")


async def xeniumranger_resegment(bundle: File, prefix: str = "xenium", args: str = "") -> File:
    """``bundle`` (a ``.tar.gz``) with its cells resegmented: ``<prefix>.tar.gz``."""
    results = await xeniumranger_resegment_cmd(bundle=bundle, prefix=prefix, args=args)
    return await xeniumranger_outs(results, prefix, "resegment")


async def xeniumranger_importsegmentation(
    bundle: File,
    transcript_assignment: File | None = None,
    viz_polygons: File | None = None,
    nuclei: File | None = None,
    cells: File | None = None,
    coordinate_transform: File | None = None,
    units: str = "",
    expansion_distance: float | None = None,
    prefix: str = "xenium",
    args: str = "",
) -> File:
    """``bundle`` (a ``.tar.gz``) with an imported segmentation: ``<prefix>.tar.gz``.

    Give either image-based results (``nuclei`` and/or ``cells``) or
    transcript-based ones (``transcript_assignment`` and ``viz_polygons``).
    """
    if (nuclei or cells) and (transcript_assignment or viz_polygons):
        raise ValueError(
            "nuclei and cells (image segmentation results) are mutually exclusive with "
            "transcript_assignment and viz_polygons (transcript assignment results)"
        )
    results = await xeniumranger_importsegmentation_cmd(
        bundle=bundle,
        transcript_assignment=transcript_assignment,
        viz_polygons=viz_polygons,
        nuclei=nuclei,
        cells=cells,
        coordinate_transform=coordinate_transform,
        units=units,
        expansion_distance="" if expansion_distance is None else f"{expansion_distance:g}",
        prefix=prefix,
        args=args,
    )
    return await xeniumranger_outs(results, prefix, "import-segmentation")
