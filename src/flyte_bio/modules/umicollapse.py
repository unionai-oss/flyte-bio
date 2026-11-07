"""umicollapse — UMI deduplication with UMICollapse.

Exposes :func:`umicollapse`, which collapses reads sharing a UMI, either in a
coordinate-sorted BAM (``mode="bam"``, by UMI and mapping position) or in a
FASTQ whose read names carry the UMI (``mode="fastq"``).

Outputs keep upstream's names: ``<prefix>.bam`` or ``<prefix>.fastq.gz``, and
the log ``<prefix>_UMICollapse.log``.

As upstream, UMICollapse's jar is run directly so its stack can be set (999
MiB; most JVMs refuse more than 1 GiB). Its heap is 80% of the task's memory
(upstream gives it 90%).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

UMICOLLAPSE_IMAGE = "quay.io/biocontainers/umicollapse:1.1.0--hdfd78af_0"
# Upstream labels it process_high and process_high_memory; UMICollapse holds a
# position's (or, for FASTQ, all) reads in memory, so large inputs need more.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="16Gi")
MODES = ("bam", "fastq")

# The BAM's index, when given, is linked next to it.
umicollapse_cmd = shell.create(
    name="umicollapse",
    image=UMICOLLAPSE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"input": File, "bai": File | None, "mode": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        IN=({{inputs.input}})
        shopt -s nullglob; BAI=({{inputs.bai}}); shopt -u nullglob
        MODE={{inputs.mode}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{IN[0]}}")
        ln -s "${{IN[0]}}" "$NAME"
        if [ ${{#BAI[@]}} -gt 0 ]; then ln -s "${{BAI[0]}}" "$(basename "${{BAI[0]}}")"; fi
        if [ "$MODE" = fastq ]; then EXT=fastq.gz; else EXT=bam; fi

        # The generated launcher allows configuring heap size, but not stack size.
        SHARE="$(dirname "$(command -v umicollapse)")/../share"
        UMICOLLAPSE_JAR=$(find "$SHARE" -maxdepth 2 -name umicollapse.jar -print -quit)
        java \
            -Xmx{JAVA_HEAP_MB}M \
            -Xss999M \
            -jar "$UMICOLLAPSE_JAR" \
            "$MODE" \
            -i "$NAME" \
            -o "$OUT/$P.$EXT" \
            $ARGS | tee "$OUT/${{P}}_UMICollapse.log"
    """,
)

env = flyte.TaskEnvironment.from_task("umicollapse", umicollapse_cmd.as_task())


@dataclass
class UmicollapseResult:
    log: File  # <prefix>_UMICollapse.log
    bam: File | None = None  # <prefix>.bam, in bam mode
    fastq: File | None = None  # <prefix>.fastq.gz, in fastq mode


async def umicollapse(
    input: File,
    mode: str,
    bai: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> UmicollapseResult:
    """Deduplicate ``input`` by UMI: a sorted BAM (``mode="bam"``) or a FASTQ (``mode="fastq"``)."""
    if mode not in MODES:
        raise ValueError(f"umicollapse mode must be one of {MODES}, got {mode!r}")
    results = await umicollapse_cmd(input=input, bai=bai, mode=mode, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"umicollapse wrote no {name}")
        return found

    result = UmicollapseResult(log=await pick(f"{prefix}_UMICollapse.log"))
    if mode == "bam":
        result.bam = await pick(f"{prefix}.bam")
    else:
        result.fastq = await pick(f"{prefix}.fastq.gz")
    return result
