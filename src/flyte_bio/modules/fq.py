"""fq — FASTQ utilities.

Exposes :func:`fq_subsample`, which samples one sample's reads (one or two
mates, kept in sync) by probability (``-p``) or record count (``-n``); one
of those must be given in ``args``. Outputs are named as upstream:
``<prefix>.fastq.gz`` (single-end) or ``<prefix>_R1/_R2.fastq.gz``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/fq:...` (Docker Hub) — use the quay.io mirror.
FQ_IMAGE = "quay.io/biocontainers/fq:0.12.0--h9ee0642_0"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# `reads_2` is `list[File]` (0 or 1 item) rather than `File | None` until
# flyteorg/flyte#8118 is deployed: copilot stages a set optional File as a
# bare path, not the per-input dir the shell glob expects.
fq_subsample_cmd = shell.create(
    name="fq_subsample",
    image=FQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": list[File], "prefix": str, "args": str},
    defaults={"reads_2": []},
    outputs={"results": Dir},
    script=r"""
        shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        if [ ${#R2[@]} -eq 0 ]; then
            fq subsample $ARGS {inputs.reads_1} --r1-dst {outputs.results}/"$P.fastq.gz"
        else
            fq subsample $ARGS {inputs.reads_1} "${R2[0]}" \
                --r1-dst {outputs.results}/"${P}_R1.fastq.gz" --r2-dst {outputs.results}/"${P}_R2.fastq.gz"
        fi
    """,
)


env = flyte.TaskEnvironment.from_task(
    "fq",
    fq_subsample_cmd.as_task(),
)


async def fq_subsample(
    reads_1: File, reads_2: File | None, args: str, prefix: str = "subsampled"
) -> tuple[File, File | None]:
    """Subsample one sample's reads; returns ``(reads_1, reads_2 or None)``."""
    if not any(flag in args.split() for flag in ("-p", "--probability", "-n", "--record-count")):
        raise ValueError("fq_subsample needs --probability (-p) or --record-count (-n) in args")
    results = await fq_subsample_cmd(
        reads_1=reads_1, reads_2=[reads_2] if reads_2 is not None else [], prefix=prefix, args=args
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"fq subsample wrote no {name}")
        return found

    if reads_2 is None:
        return await pick(f"{prefix}.fastq.gz"), None
    return await pick(f"{prefix}_R1.fastq.gz"), await pick(f"{prefix}_R2.fastq.gz")
