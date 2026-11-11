"""datavzrd — Datavzrd, interactive HTML reports of tabular data.

Exposes :func:`datavzrd`, which renders the tables a YAML config describes
(CSV, TSV …) into a self-contained interactive HTML report: upstream's
``<prefix>`` directory (``index.html``, ``static/`` and one directory per
view). The config names its tables by file name: they're linked into the
working directory beside it, as upstream stages them.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

DATAVZRD_IMAGE = "community.wave.seqera.io/library/datavzrd:2.75.3--58fb85a496607c2e"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

datavzrd_cmd = shell.create(
    name="datavzrd",
    image=DATAVZRD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"config": File | None, "files": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        shopt -s nullglob; CONFIG=({inputs.config}); FILES=({inputs.files}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        if [ ${#CONFIG[@]} -eq 0 ]; then
            echo "datavzrd needs a config file" >&2
            exit 1
        fi
        W=$(mktemp -d); cd "$W"
        for f in "${CONFIG[0]}" "${FILES[@]}"; do ln -s "$f" "$(basename "$f")"; done
        mkdir "$P"
        # datavzrd is a Rust tool: a panic's exit code could be retried forever, so any failure exits 1.
        datavzrd $ARGS "$(basename "${CONFIG[0]}")" --output "$P" || exit 1
        mv "$P" {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("datavzrd", datavzrd_cmd.as_task())


async def datavzrd(config: File | None, files: list[File], prefix: str, args: str = "") -> Dir:
    """The HTML report of ``files`` that ``config`` describes: the ``<prefix>`` Dir.

    Fails without a ``config``, as upstream.
    """
    results = await datavzrd_cmd(config=config, files=files, prefix=prefix, args=args)
    return Dir.from_existing_remote(results.path.rstrip("/") + f"/{prefix}")
