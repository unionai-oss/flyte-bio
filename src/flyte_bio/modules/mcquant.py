"""mcquant — single-cell quantification of a multiplexed image over a segmentation mask.

Exposes :func:`mcquant`, which measures, for every cell (label) in a
segmentation mask, each channel's intensity (mean by default; ``args`` such as
``--intensity_props intensity_sum`` add more) and the cell's morphology
(centroid, area, axis lengths, eccentricity, solidity, extent, orientation).
Channels are named from the markers CSV's ``marker_name`` column. The table is
written as ``<image stem>_<mask stem>.csv`` (stems drop ``.ome.tif`` / ``.tif``),
e.g. ``cycif_tonsil_registered_cell.csv``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MCQUANT_IMAGE = "quay.io/nf-core/quantification:1.5.4"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

mcquant_cmd = shell.create(
    name="mcquant",
    image=MCQUANT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "mask": File, "markers": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGE=({{inputs.image}})
        MASK=({{inputs.mask}})
        MARKERS=({{inputs.markers}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib"
        python /app/CommandSingleCellExtraction.py \
            --masks "${{MASK[0]}}" \
            --image "${{IMAGE[0]}}" \
            --channel_names "${{MARKERS[0]}}" \
            --output . \
            $ARGS || exit 1
        shopt -s nullglob
        TABLES=(./*.csv)
        [ ${{#TABLES[@]}} -gt 0 ] || {{ echo "mcquant wrote no CSV" >&2; exit 1; }}
        mv "${{TABLES[@]}}" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("mcquant", mcquant_cmd.as_task())


async def mcquant(image: File, mask: File, markers: File, args: str = "") -> File:
    """The cell-by-feature table of ``image`` over ``mask``: ``<image stem>_<mask stem>.csv``."""
    results = await mcquant_cmd(image=image, mask=mask, markers=markers, args=args)
    tables = [f for f in await results.list_files() if f.name.endswith(".csv")]
    if len(tables) != 1:
        raise FileNotFoundError(f"mcquant wrote {len(tables)} CSVs, expected 1: {[t.name for t in tables]}")
    return tables[0]
