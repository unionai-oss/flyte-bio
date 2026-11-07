"""pear — PEAR, merging of overlapping paired-end reads.

Exposes :func:`pear`, which merges one sample's read pairs and returns, under
upstream's names:

- ``<prefix>.assembled.fastq.gz``: the merged reads;
- ``<prefix>.unassembled.forward.fastq.gz`` and
  ``<prefix>.unassembled.reverse.fastq.gz``: the pairs that didn't merge;
- ``<prefix>.discarded.fastq.gz``: the reads PEAR discarded.

As upstream, the reads are decompressed first (into the working directory,
under their names without ``.gz``), and PEAR runs on the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PEAR_IMAGE = "quay.io/biocontainers/pear:0.9.6--h67092d7_8"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

pear_cmd = shell.create(
    name="pear",
    image=PEAR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        R2=({{inputs.reads_2}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        F=$(basename "${{R1[0]}}" .gz)
        R=$(basename "${{R2[0]}}" .gz)
        gzip -cd "${{R1[0]}}" > "$F"
        gzip -cd "${{R2[0]}}" > "$R"
        pear -f "$F" -r "$R" -o "$P" -j $TASK_CPUS $ARGS || exit 1
        for KIND in assembled unassembled.forward unassembled.reverse discarded; do
            gzip -f "$P.$KIND.fastq"
            mv "$P.$KIND.fastq.gz" "$OUT"/
        done
    """,
)

env = flyte.TaskEnvironment.from_task("pear", pear_cmd.as_task())


@dataclass
class PearResult:
    assembled: File  # <prefix>.assembled.fastq.gz
    unassembled: list[File]  # <prefix>.unassembled.{forward,reverse}.fastq.gz
    discarded: File  # <prefix>.discarded.fastq.gz
    results: Dir


async def pear(reads_1: File, reads_2: File, prefix: str = "sample", args: str = "") -> PearResult:
    """Merge the read pairs in ``reads_1`` and ``reads_2``."""
    results = await pear_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)

    async def pick(kind: str) -> File:
        name = f"{prefix}.{kind}.fastq.gz"
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"pear wrote no {name}")
        return found

    return PearResult(
        assembled=await pick("assembled"),
        unassembled=[await pick("unassembled.forward"), await pick("unassembled.reverse")],
        discarded=await pick("discarded"),
        results=results,
    )
