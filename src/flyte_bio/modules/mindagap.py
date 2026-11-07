"""mindagap — fill the grid lines of tiled spatial-transcriptomics panoramas, and mark duplicated spots.

Exposes:

- :func:`mindagap_mindagap`: fill the empty grid lines between the tiles of a
  panorama image (TIFF or PNG) with neighbour-weighted values (Mindagap's
  ``mindagap.py``), writing ``<image name>_gridfilled.<ext>``.
- :func:`mindagap_duplicatefinder`: mark the transcript spots duplicated in the
  overlap of neighbouring tiles in an XYZ spot table (Mindagap's
  ``duplicate_finder.py``), writing ``<table name>_markedDups.txt``.

Both take Mindagap's options in ``args`` (e.g. ``"3 40 --Xtilesize 2144"`` or
``"90"``, as upstream's tests). Mindagap writes its output next to its input,
so the input is linked into a working directory first.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MINDAGAP_IMAGE = "quay.io/biocontainers/mindagap:0.0.2--pyhdfd78af_1"
# Upstream labels mindagap process_low and duplicatefinder process_single; large panoramas need more memory.
MINDAGAP_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")
DUPLICATEFINDER_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

mindagap_mindagap_cmd = shell.create(
    name="mindagap_mindagap",
    image=MINDAGAP_IMAGE,
    resources=MINDAGAP_RESOURCES,
    inputs={"panorama": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PANORAMA=({inputs.panorama})
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib"
        NAME=$(basename "${PANORAMA[0]}")
        ln -s "${PANORAMA[0]}" "$NAME"
        mindagap.py "$NAME" $ARGS || exit 1
        shopt -s nullglob
        FILLED=(*_gridfilled.tif *_gridfilled.tiff *_gridfilled.png)
        if [ ${#FILLED[@]} -eq 0 ]; then echo "mindagap wrote no gridfilled image" >&2; exit 1; fi
        mv "${FILLED[@]}" "$OUT"/
    """,
)

mindagap_duplicatefinder_cmd = shell.create(
    name="mindagap_duplicatefinder",
    image=MINDAGAP_IMAGE,
    resources=DUPLICATEFINDER_RESOURCES,
    inputs={"spot_table": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SPOTS=({inputs.spot_table})
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib"
        NAME=$(basename "${SPOTS[0]}")
        ln -s "${SPOTS[0]}" "$NAME"
        duplicate_finder.py "$NAME" $ARGS || exit 1
        mv "${NAME%.*}_markedDups.txt" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "mindagap", mindagap_mindagap_cmd.as_task(), mindagap_duplicatefinder_cmd.as_task()
)


async def mindagap_mindagap(panorama: File, args: str = "") -> File:
    """``panorama`` with its tile grid lines filled: ``<panorama name>_gridfilled.<ext>``."""
    results = await mindagap_mindagap_cmd(panorama=panorama, args=args)
    files = await results.list_files()
    assert len(files) == 1, f"mindagap wrote {len(files)} gridfilled images, expected 1"
    return files[0]


async def mindagap_duplicatefinder(spot_table: File, args: str = "") -> File:
    """``spot_table`` with the spots duplicated across tile overlaps marked: ``<table name>_markedDups.txt``."""
    results = await mindagap_duplicatefinder_cmd(spot_table=spot_table, args=args)
    files = await results.list_files()
    assert len(files) == 1, f"duplicate_finder wrote {len(files)} tables, expected 1"
    return files[0]
