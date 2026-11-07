"""fastqutils — fastq_utils FASTQ validation.

Exposes :func:`fastqutils_info`, which validates one sample's FASTQ (or a
pair of mates) with ``fastq_info``. A malformed file (truncated, sequence and
quality lengths differing, unpaired mates) makes the task fail with
fastq_info's error. On success, as upstream, it writes ``<prefix>.txt``
naming the files that were checked.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

FASTQ_UTILS_IMAGE = "quay.io/biocontainers/fastq_utils:0.25.3--ha9dfd29_0"
# Upstream labels fastq_info process_medium (6 CPUs, 36 GB); it is
# single-threaded, and its read-name index stays small for test-sized files.
INFO_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# The reads are linked into a working directory under their own names, so the
# report names them as upstream's does. On failure fastq_info's error is also
# written to the container's termination message, so it reaches the task's
# error (the captured stderr isn't uploaded from a failed task).
fastqutils_info_cmd = shell.create(
    name="fastqutils_info",
    image=FASTQ_UTILS_IMAGE,
    resources=INFO_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        R1=({inputs.reads_1})
        shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        READS=()
        for F in "${R1[@]}" "${R2[@]}"; do ln -s "$F" "$(basename "$F")"; READS+=("$(basename "$F")"); done
        if ! fastq_info $ARGS "${READS[@]}" 2> fastq_info.err; then
            cat fastq_info.err >&2
            { tail -c 2048 fastq_info.err > /dev/termination-log; } 2> /dev/null
            exit 1
        fi
        cat fastq_info.err >&2
        echo "fastq_utils fastq_info ran and found no issues with ${READS[*]}" > {outputs.results}/"$P".txt
    """,
)

env = flyte.TaskEnvironment.from_task("fastqutils", fastqutils_info_cmd.as_task())


async def fastqutils_info(
    reads_1: File, reads_2: File | None = None, prefix: str = "fastqutils", args: str = ""
) -> File:
    """Validate ``reads_1`` (and its mate ``reads_2``); fails on a malformed file, else returns ``<prefix>.txt``."""
    results = await fastqutils_info_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)
    report = await results.get_file(f"{prefix}.txt")
    if report is None:
        raise FileNotFoundError(f"fastq_info wrote no {prefix}.txt")
    return report
