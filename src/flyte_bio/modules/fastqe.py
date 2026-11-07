"""fastqe — FASTQ quality summaries as emoji.

Exposes :func:`fastqe`, which summarizes the base qualities of one or more
FASTQ files into ``<prefix>.tsv``: per file, the mean (and with ``--min`` /
``--max`` in ``args``, minimum and maximum) quality at each position, as
emoji.

The TSV names each input file as given on the command line, so, as upstream,
the files are passed by their bare names.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

FASTQE_IMAGE = "quay.io/biocontainers/fastqe:0.5.2--pyhdfd78af_0"

# Upstream's process_single: one CPU.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# List inputs are staged in file-name order.
fastqe_cmd = shell.create(
    name="fastqe",
    image=FASTQE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        READS=({inputs.reads})
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for r in "${READS[@]}"; do ln -s "$r" .; NAMES+=("$(basename "$r")"); done
        fastqe $ARGS "${NAMES[@]}" --output {outputs.results}/{inputs.prefix}.tsv || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("fastqe", fastqe_cmd.as_task())


async def fastqe(reads: list[File], prefix: str, args: str = "") -> File:
    """fastqe's emoji quality summary of ``reads``: ``<prefix>.tsv``."""
    results = await fastqe_cmd(reads=reads, prefix=prefix, args=args)
    tsv = await results.get_file(f"{prefix}.tsv")
    if tsv is None:
        raise FileNotFoundError(f"fastqe wrote no {prefix}.tsv")
    return tsv
