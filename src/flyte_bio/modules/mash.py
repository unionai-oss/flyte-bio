"""mash — fast genome and metagenome distances with MinHash sketches.

Exposes:

- :func:`mash_sketch` — sketch sequences (``mash sketch``): one sketch per
  input file, all in ``<prefix>.msh``, with Mash's report of the sketching
  (its standard error) as ``<prefix>.mash_stats``.
- :func:`mash_dist` — the distances between the sequences or sketches of a
  query and a reference (``mash dist``), as ``<prefix>.txt``.
- :func:`mash_screen` — which of a sketch's sequences are contained in a
  query's reads or sequences (``mash screen``), as ``<prefix>.screen``.

Mash records and prints its inputs' paths, so they're staged in a working
directory under their own names and passed by name. Several inputs are taken
in their names' alphabetical order. ``-p`` is the task's CPUs; options go in
``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MASH_IMAGE = "quay.io/biocontainers/mash:2.3--he348c14_1"
# Upstream: sketch and screen process_medium, dist process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
DIST_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# Links each of the given files into the working directory under its own name.
STAGE_BY_NAME = r"""
        stage_by_name() {
            local F
            for F in "$@"; do ln -sf "$F" "$(basename "$F")"; done
        }
""".strip()

mash_sketch_cmd = shell.create(
    name="mash_sketch",
    image=MASH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {STAGE_BY_NAME}
        IN=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        stage_by_name "${{IN[@]}}"
        NAMES=()
        for F in "${{IN[@]}}"; do NAMES+=("$(basename "$F")"); done
        mash \
            sketch \
            $ARGS \
            "${{NAMES[@]}}" \
            -p $TASK_CPUS \
            -o "$P" \
            2> "$P.mash_stats" \
            || {{ cat "$P.mash_stats" >&2; exit 1; }}
        cat "$P.mash_stats" >&2
        mv "$P.msh" "$P.mash_stats" "$OUT"/
    """,
)

mash_dist_cmd = shell.create(
    name="mash_dist",
    image=MASH_IMAGE,
    resources=DIST_RESOURCES,
    inputs={"query": File, "reference": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {STAGE_BY_NAME}
        QUERY=({{inputs.query}}); REF=({{inputs.reference}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        stage_by_name "${{QUERY[0]}}" "${{REF[0]}}"
        mash \
            dist \
            -p $TASK_CPUS \
            $ARGS \
            "$(basename "${{REF[0]}}")" \
            "$(basename "${{QUERY[0]}}")" > "$OUT/$P.txt" \
            || exit 1
    """,
)

mash_screen_cmd = shell.create(
    name="mash_screen",
    image=MASH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"query": list[File], "sketch": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {STAGE_BY_NAME}
        IN=({{inputs.query}}); SKETCH=({{inputs.sketch}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        stage_by_name "${{IN[@]}}" "${{SKETCH[0]}}"
        NAMES=()
        for F in "${{IN[@]}}"; do NAMES+=("$(basename "$F")"); done
        mash \
            screen \
            $ARGS \
            -p $TASK_CPUS \
            "$(basename "${{SKETCH[0]}}")" \
            "${{NAMES[@]}}" \
            > "$OUT/$P.screen" \
            || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "mash", mash_sketch_cmd.as_task(), mash_dist_cmd.as_task(), mash_screen_cmd.as_task()
)


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def mash_sketch(reads: File | list[File], prefix: str = "sample", args: str = "") -> tuple[File, File]:
    """Sketch each of ``reads``: the sketches ``<prefix>.msh`` and the report ``<prefix>.mash_stats``."""
    inputs = reads if isinstance(reads, list) else [reads]
    results = await mash_sketch_cmd(reads=inputs, prefix=prefix, args=args)
    return (
        await pick_file(results, f"{prefix}.msh", "mash sketch"),
        await pick_file(results, f"{prefix}.mash_stats", "mash sketch"),
    )


async def mash_dist(query: File, reference: File, prefix: str = "sample", args: str = "") -> File:
    """The Mash distances between ``query`` and ``reference`` (sequences or sketches): ``<prefix>.txt``."""
    results = await mash_dist_cmd(query=query, reference=reference, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.txt", "mash dist")


async def mash_screen(query: File | list[File], sketch: File, prefix: str = "sample", args: str = "") -> File:
    """Which of ``sketch``'s sequences ``query``'s reads or sequences contain: ``<prefix>.screen``."""
    inputs = query if isinstance(query, list) else [query]
    results = await mash_screen_cmd(query=inputs, sketch=sketch, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.screen", "mash screen")
