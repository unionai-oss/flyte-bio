"""fusioncatcher — gene-fusion detection from RNA-seq reads with FusionCatcher.

Exposes :func:`fusioncatcher`, which finds fusion genes in one sample's reads
against a FusionCatcher database, writing, under upstream's names:

- ``<prefix>.fusion-genes.txt`` (``final-list_candidate-fusion-genes.txt``);
- ``<prefix>.summary.txt`` (``summary_candidate_fusions.txt``);
- ``<prefix>.log`` (``fusioncatcher.log``).

The database (``human_v102`` and the like, about 50 GB unpacked) is a
directory far too big for Flyte's input staging, so it travels as a
``.tar.gz`` of the directory (as FusionCatcher distributes it), extracted in
the task: the archive's single top-level directory, if it has one, is the
database. Unpacking needs that much scratch disk. The Java heap (``--Xmx``) is
80% of the task's memory, as upstream; a human database needs more memory than
the default (override ``resources``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

FUSIONCATCHER_IMAGE = "quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_5"
# Upstream's process_medium, capped.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="24Gi")

# As upstream: the reads are passed comma-separated (R1 first), and the
# outputs renamed to the prefix.
fusioncatcher_cmd = shell.create(
    name="fusioncatcher",
    image=FUSIONCATCHER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "reference": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}}); REFERENCE=({{inputs.reference}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        DATA_ROOT=$(mktemp -d)
        tar -xzf "${{REFERENCE[0]}}" -C "$DATA_ROOT" || exit 1
        DATA_ENTRIES=("$DATA_ROOT"/*)
        DATA="$DATA_ROOT"
        if [ ${{#DATA_ENTRIES[@]}} -eq 1 ] && [ -d "${{DATA_ENTRIES[0]}}" ]; then DATA="${{DATA_ENTRIES[0]}}"; fi
        INPUT="${{R1[0]}}"
        if [ ${{#R2[@]}} -gt 0 ]; then INPUT="$INPUT,${{R2[0]}}"; fi
        W=$(mktemp -d); cd "$W"
        fusioncatcher \
            --input="$INPUT" \
            --output=. \
            --data="$DATA" \
            --threads=$TASK_CPUS \
            --Xmx={JAVA_HEAP_MB}M \
            $ARGS || exit 1
        mv final-list_candidate-fusion-genes.txt "$OUT/$P.fusion-genes.txt"
        mv summary_candidate_fusions.txt "$OUT/$P.summary.txt"
        mv fusioncatcher.log "$OUT/$P.log"
    """,
)

env = flyte.TaskEnvironment.from_task("fusioncatcher", fusioncatcher_cmd.as_task())


@dataclass
class FusionCatcherResult:
    fusions: File | None  # <prefix>.fusion-genes.txt
    summary: File | None  # <prefix>.summary.txt
    log: File  # <prefix>.log


async def fusioncatcher(
    reads_1: File, reference: File, prefix: str, reads_2: File | None = None, args: str = ""
) -> FusionCatcherResult:
    """Fusion genes in one sample's reads (``reads_2`` None for single-end) against the database archive.

    ``reference`` is a FusionCatcher database directory as a ``.tar.gz``.
    """
    results = await fusioncatcher_cmd(reads_1=reads_1, reads_2=reads_2, reference=reference, prefix=prefix, args=args)
    log = await results.get_file(f"{prefix}.log")
    if log is None:
        raise FileNotFoundError(f"fusioncatcher wrote no {prefix}.log")
    return FusionCatcherResult(
        fusions=await results.get_file(f"{prefix}.fusion-genes.txt"),
        summary=await results.get_file(f"{prefix}.summary.txt"),
        log=log,
    )
