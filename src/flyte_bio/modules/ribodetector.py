"""ribodetector — remove rRNA reads from FASTQ with a deep-learning classifier.

Exposes :func:`ribodetector`, which runs the CPU build (``ribodetector_cpu``)
on one sample's reads and keeps the non-rRNA reads, as upstream without an
accelerator: ``<prefix>.nonrna.fastq.gz`` (single-end) or
``<prefix>.nonrna.1.fastq.gz`` / ``<prefix>.nonrna.2.fastq.gz`` (paired-end),
plus the log ``<prefix>.log``. ``length`` is the reads' mean length (``-l``).
Upstream's GPU variant (another image, ``ribodetector -m``) isn't ported.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RIBODETECTOR_IMAGE = "community.wave.seqera.io/library/ribodetector:0.3.3--ad3d7071e408b502"
# Upstream labels it process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

ribodetector_cmd = shell.create(
    name="ribodetector",
    image=RIBODETECTOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "length": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        READS=("$(basename "${{R1[0]}}")"); ln -s "${{R1[0]}}" "${{READS[0]}}"
        if [ ${{#R2[@]}} -gt 0 ]; then
            READS+=("$(basename "${{R2[0]}}")"); ln -s "${{R2[0]}}" "${{READS[1]}}"
            OUTPUTS=("$P".nonrna.1.fastq.gz "$P".nonrna.2.fastq.gz)
        else
            OUTPUTS=("$P".nonrna.fastq.gz)
        fi
        ribodetector_cpu \
            -i "${{READS[@]}}" \
            -o "${{OUTPUTS[@]}}" \
            -l {{inputs.length}} \
            -t $TASK_CPUS \
            --log "$P".log \
            $ARGS \
            || exit 1
        mv "${{OUTPUTS[@]}}" "$P".log "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("ribodetector", ribodetector_cmd.as_task())


@dataclass
class RibodetectorResult:
    reads: list[File]  # <prefix>.nonrna.fastq.gz, or <prefix>.nonrna.{1,2}.fastq.gz
    log: File  # <prefix>.log


async def ribodetector(
    reads_1: File, length: int, reads_2: File | None = None, prefix: str = "sample", args: str = ""
) -> RibodetectorResult:
    """Keep the non-rRNA reads of one sample (``reads_2`` None for single-end) of mean length ``length``."""
    results = await ribodetector_cmd(reads_1=reads_1, reads_2=reads_2, length=length, prefix=prefix, args=args)
    names = [f"{prefix}.nonrna.fastq.gz"] if reads_2 is None else [f"{prefix}.nonrna.{i}.fastq.gz" for i in (1, 2)]

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"ribodetector wrote no {name}")
        return found

    return RibodetectorResult(reads=[await pick(n) for n in names], log=await pick(f"{prefix}.log"))
