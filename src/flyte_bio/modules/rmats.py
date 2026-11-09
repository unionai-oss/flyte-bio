"""rmats — rMATS alternative splicing, its ``prep`` step.

Exposes :func:`rmats_prep`, which runs ``rmats.py --task prep`` on one BAM
against a GTF and writes ``<prefix>.rmats``, the per-BAM splice graph and
junction counts that rMATS' ``post`` step combines, and
``<prefix>_read_outcomes_by_bam.txt``, the counts of reads used and of each
reason for skipping one.

As upstream, the BAM is listed by its name in ``<prefix>.prep.b1.txt`` and
rMATS runs on the task's CPUs. ``args`` takes rMATS' options, e.g.
``--variable-read-length --allow-clipping``, ``-t single`` or
``--libType fr-firststrand``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RMATS_IMAGE = "community.wave.seqera.io/library/rmats:4.3.0--85c9d01519f77bc5"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

rmats_prep_cmd = shell.create(
    name="rmats_prep",
    image=RMATS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "gtf": File, "read_length": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); GTF=({{inputs.gtf}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}"); GTF_NAME=$(basename "${{GTF[0]}}")
        ln -s "${{BAM[0]}}" "$BAM_NAME"; ln -s "${{GTF[0]}}" "$GTF_NAME"
        echo "$BAM_NAME" > "$P.prep.b1.txt"
        rmats.py --task prep $ARGS --nthread $TASK_CPUS --b1 "$P.prep.b1.txt" --gtf "$GTF_NAME" \
            --readLength {{inputs.read_length}} --tmp "${{P}}_rmats_tmp" --od "${{P}}_rmats_prep" || exit 1
        cp "${{P}}"_rmats_tmp/*.txt "$OUT/${{P}}_read_outcomes_by_bam.txt"
        cp "${{P}}"_rmats_tmp/*.rmats "$OUT/$P.rmats"
    """,
)

env = flyte.TaskEnvironment.from_task("rmats", rmats_prep_cmd.as_task())


@dataclass
class RmatsPrepResult:
    rmats: File  # <prefix>.rmats
    read_outcomes: File  # <prefix>_read_outcomes_by_bam.txt


async def rmats_prep(bam: File, gtf: File, read_length: int, prefix: str, args: str = "") -> RmatsPrepResult:
    """rMATS' prep step for ``bam`` against ``gtf``, with reads of ``read_length``."""
    results = await rmats_prep_cmd(bam=bam, gtf=gtf, read_length=read_length, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"rmats_prep wrote no {name}")
        return found

    return RmatsPrepResult(
        rmats=await pick(f"{prefix}.rmats"), read_outcomes=await pick(f"{prefix}_read_outcomes_by_bam.txt")
    )
