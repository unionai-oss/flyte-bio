"""find — batch file operations driven by ``find``.

Exposes :func:`find_unpigz`, which decompresses gzipped files with
``unpigz``: every file is linked into ``gzipped/``, found with ``find`` and
decompressed to ``ungzipped/<name without .gz>``, as upstream. Every input
must end in ``.gz``. ``--processes`` is the task's CPUs; options go in
``args``. The decompressed files come back in one Dir, so pass at most a few
hundred files per call.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FIND_UNPIGZ_IMAGE = "community.wave.seqera.io/library/findutils_pigz:c4dd5edc44402661"
# Upstream labels it process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="4Gi")

NOT_GZIPPED = "All files provided to this module must be gzipped (and have the .gz extension)."

find_unpigz_cmd = shell.create(
    name="find_unpigz",
    image=FIND_UNPIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"files": list[File], "args": str},
    defaults={"args": ""},
    outputs={"ungzipped": Dir},
    script=rf"""
        {TASK_LIMITS}
        FILES=({{inputs.files}})
        ARGS={{inputs.args}}
        OUT={{outputs.ungzipped}}
        W=$(mktemp -d); cd "$W"
        mkdir gzipped
        for F in "${{FILES[@]}}"; do ln -s "$F" gzipped/"$(basename "$F")"; done
        while IFS= read -r -d $'\0' file; do
            unpigz \
                $ARGS \
                -cd \
                --processes $TASK_CPUS \
                "$file" \
                > "$OUT"/"$(basename "$file" .gz)" \
                || exit 1
        done < <( find gzipped/ -name '*.gz' -print0 )
    """,
)

env = flyte.TaskEnvironment.from_task("find", find_unpigz_cmd.as_task())


async def find_unpigz(files: File | list[File], args: str = "") -> list[File]:
    """Decompress gzipped ``files``: each named like its input without ``.gz``, in input order."""
    gzipped = files if isinstance(files, list) else [files]
    if not gzipped or any(not f.name.endswith(".gz") for f in gzipped):
        raise ValueError(NOT_GZIPPED)
    results = await find_unpigz_cmd(files=gzipped, args=args)
    out = []
    for f in gzipped:
        name = f.name.removesuffix(".gz")
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"find_unpigz wrote no {name}")
        out.append(found)
    return out
