"""mat2json — convert MATLAB ``.mat`` files to JSON or CSV.

Exposes :func:`mat2json`, which converts one ``.mat`` file with the mat2json
tool and returns what it wrote, ``<mat name>.json`` and/or ``<mat name>.csv``,
sorted by name. Upstream files them under ``<process>/<prefix>/``; here they
are returned as files.

The tool is a compiled MATLAB application. Its image (MIT, like the tool)
ships the MATLAB Runtime R2023a, which MathWorks distributes under its own
MATLAB Runtime licence: free to use and redistribute, but not open source.
"""

import os

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MAT2JSON_IMAGE = "quay.io/nf-core/mat2json:1.0.0"
# Upstream labels it process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

mat2json_cmd = shell.create(
    name="mat2json",
    image=MAT2JSON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"matfile": File},
    outputs={"results": Dir},
    script=r"""
        MAT=({inputs.matfile})
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        # The MATLAB Runtime unpacks the application into its cache.
        export MCR_CACHE_ROOT="$W/.mcr" HOME="$W"
        mkdir -p "$MCR_CACHE_ROOT"
        # mat2json writes its outputs next to the .mat file.
        ln -s "${MAT[0]}" .
        mat2json "$(basename "${MAT[0]}")" || exit 1
        shopt -s nullglob
        CONVERTED=(*.json *.csv)
        if [ ${#CONVERTED[@]} -eq 0 ]; then echo "mat2json wrote no .json or .csv" >&2; exit 1; fi
        mv "${CONVERTED[@]}" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("mat2json", mat2json_cmd.as_task())


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def mat2json(matfile: File) -> list[File]:
    """Convert ``matfile`` to JSON or CSV: ``<name>.json`` and/or ``<name>.csv``, sorted by name."""
    results = await mat2json_cmd(matfile=matfile)
    return sorted(await results.list_files(), key=file_name)
