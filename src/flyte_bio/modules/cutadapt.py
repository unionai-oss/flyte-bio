"""cutadapt — adapter, primer and quality trimming of FASTQ reads.

Exposes :func:`cutadapt`, which trims one sample's reads (single-end, or
paired-end with ``reads_2``) and returns the trimmed reads and cutadapt's log.
The adapters and trimming options come in ``args`` (upstream's ``ext.args``),
for example ``-a AGATCGGAAGAGC -q 25``.

Outputs keep upstream's names:

- trimmed reads: ``<prefix>.trim.fastq.gz`` (single-end) or
  ``<prefix>_1.trim.fastq.gz`` / ``<prefix>_2.trim.fastq.gz`` (paired-end);
- the report cutadapt prints to stdout: ``<prefix>.cutadapt.log``.

As upstream, ``--cores`` is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CUTADAPT_IMAGE = "community.wave.seqera.io/library/cutadapt:5.2--efafca34ac3f6c8a"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="4Gi")

cutadapt_cmd = shell.create(
    name="cutadapt",
    image=CUTADAPT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        W=$(mktemp -d); cd "$W"
        if [ ${{#R2[@]}} -eq 0 ]; then
            cutadapt --cores $TASK_CPUS $ARGS -o "$P.trim.fastq.gz" "${{R1[0]}}" > "$P.cutadapt.log"
        else
            cutadapt --cores $TASK_CPUS $ARGS -o "${{P}}_1.trim.fastq.gz" -p "${{P}}_2.trim.fastq.gz" \
                "${{R1[0]}}" "${{R2[0]}}" > "$P.cutadapt.log"
        fi
        cat "$P.cutadapt.log" >&2
        set +f
        mv ./*.trim.fastq.gz "$P.cutadapt.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("cutadapt", cutadapt_cmd.as_task())


@dataclass
class CutadaptResult:
    reads: list[File]  # the trimmed reads: one file, or mates 1 and 2
    log: File  # <prefix>.cutadapt.log
    results: Dir  # everything cutadapt wrote


async def cutadapt(
    reads_1: File, reads_2: File | None = None, prefix: str = "sample", args: str = ""
) -> CutadaptResult:
    """Trim one sample's reads; ``reads_2`` is None for single-end reads."""
    results = await cutadapt_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"cutadapt wrote no {name}")
        return found

    mates = ["_1", "_2"] if reads_2 is not None else [""]
    return CutadaptResult(
        reads=[await pick(f"{prefix}{mate}.trim.fastq.gz") for mate in mates],
        log=await pick(f"{prefix}.cutadapt.log"),
        results=results,
    )
