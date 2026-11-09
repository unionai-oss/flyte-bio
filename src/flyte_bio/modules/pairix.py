"""pairix — index a block-compressed, sorted contact pairs file.

Exposes :func:`pairix`, which indexes a bgzipped pairs file (sorted by
chromosome pair, then position) and returns the pairs with their index,
``<pairs name>.px2``, as upstream. ``args`` is upstream's ``ext.args`` (e.g.
``-p pairs``, or column options for other formats) and is split on whitespace.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PAIRIX_IMAGE = "quay.io/biocontainers/pairix:0.3.7--py36h30a8e3e_3"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# pairix writes the index next to its input, so the input is linked into a
# working directory (inputs are read-only).
pairix_cmd = shell.create(
    name="pairix",
    image=PAIRIX_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"pairs": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PAIRS=({inputs.pairs})
        ARGS={inputs.args}
        NAME=$(basename "${PAIRS[0]}")
        W=$(mktemp -d); cd "$W"
        ln -s "${PAIRS[0]}" "$NAME"
        pairix $ARGS "$NAME" || exit 1
        mv "$NAME.px2" {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("pairix", pairix_cmd.as_task())


@dataclass
class PairixResult:
    pairs: File  # the input pairs, unchanged
    index: File  # <pairs name>.px2
    results: Dir


async def pairix(pairs: File, args: str = "") -> PairixResult:
    """Index ``pairs``: returns them with ``<pairs name>.px2``."""
    results = await pairix_cmd(pairs=pairs, args=args)
    indexes = [file for file in await results.list_files() if file.path.endswith(".px2")]
    if len(indexes) != 1:
        raise FileNotFoundError(f"pairix wrote {len(indexes)} .px2 indexes, expected 1")
    return PairixResult(pairs=pairs, index=indexes[0], results=results)
