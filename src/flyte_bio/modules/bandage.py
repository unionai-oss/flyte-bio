"""bandage — draws an assembly graph as an image.

Exposes :func:`bandage_image`, which renders a GFA assembly graph (plain or
gzipped, which is decompressed first, as upstream) to ``<prefix>.png`` and
``<prefix>.svg``. ``args`` carries Bandage's image options (e.g.
``--height 1000``). Bandage is a Qt program; it draws offscreen.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BANDAGE_IMAGE = "community.wave.seqera.io/library/bandage:0.9.0--4f0567049a14ea6d"
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

bandage_image_cmd = shell.create(
    name="bandage_image",
    image=BANDAGE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gfa": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GFA=({inputs.gfa})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        export QT_QPA_PLATFORM=offscreen
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${GFA[0]}")
        if [[ "$NAME" == *.gz ]]; then
            NAME=${NAME%.gz}
            zcat "${GFA[0]}" > "$NAME"
        else
            ln -s "${GFA[0]}" "$NAME"
        fi
        Bandage image "$NAME" "$P.png" $ARGS || exit 1
        Bandage image "$NAME" "$P.svg" $ARGS || exit 1
        mv "$P.png" "$P.svg" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("bandage", bandage_image_cmd.as_task())


@dataclass
class BandageImageResult:
    png: File  # <prefix>.png
    svg: File  # <prefix>.svg


async def bandage_image(gfa: File, prefix: str = "sample", args: str = "") -> BandageImageResult:
    """Images of the assembly graph ``gfa`` (``.gfa`` or ``.gfa.gz``)."""
    results = await bandage_image_cmd(gfa=gfa, prefix=prefix, args=args)
    images = []
    for name in (f"{prefix}.png", f"{prefix}.svg"):
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"Bandage wrote no {name}")
        images.append(found)
    return BandageImageResult(png=images[0], svg=images[1])
