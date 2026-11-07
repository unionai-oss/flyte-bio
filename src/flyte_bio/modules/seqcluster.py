"""seqcluster — small-RNA read processing.

Exposes :func:`seqcluster_collapse`, ``seqcluster collapse``: collapse
identical reads of one FASTQ (``args`` sets e.g. ``-m 1 --min_size 15``), as
upstream gzipped to ``<prefix>.fastq.gz``. A prefix that would overwrite the
input is an error.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SEQCLUSTER_IMAGE = "quay.io/biocontainers/seqcluster:1.2.9--pyh5e36f6f_0"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

seqcluster_collapse_cmd = shell.create(
    name="seqcluster_collapse",
    image=SEQCLUSTER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastq": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FQ=({inputs.fastq})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FQ[0]}"); ln -s "${FQ[0]}" "$NAME"

        seqcluster collapse $ARGS -f "$NAME" -o collapsed
        gzip collapsed/*_trimmed.fastq
        mv collapsed/*_trimmed.fastq.gz "$OUT/$P.fastq.gz"
    """,
)

env = flyte.TaskEnvironment.from_task("seqcluster", seqcluster_collapse_cmd.as_task())


async def seqcluster_collapse(fastq: File, prefix: str = "sample", args: str = "") -> File:
    """Collapse identical reads of ``fastq``: ``<prefix>.fastq.gz``."""
    if fastq.path.rstrip("/").rsplit("/", 1)[-1] == f"{prefix}.fastq.gz":
        raise ValueError("Input and output names are the same, set prefix in module configuration to disambiguate!")
    results = await seqcluster_collapse_cmd(fastq=fastq, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.fastq.gz")
    if found is None:
        raise FileNotFoundError(f"seqcluster collapse wrote no {prefix}.fastq.gz")
    return found
