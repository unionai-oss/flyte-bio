"""fastqc — read quality-control reports.

Exposes :data:`fastqc`, which writes FastQC's HTML report and data zip for
one sample's reads (one or two mates) into an output Dir, and
:data:`fastqc_reads`, the same for any number of FASTQs (e.g. every lane of
a sample).

As upstream, the reads are first symlinked to ``<prefix>.gz`` (single-end)
or ``<prefix>_1.gz`` / ``<prefix>_2.gz`` (paired) so the reports are named
after the sample, e.g. ``<prefix>_1_fastqc.html``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/fastqc:...` (Docker Hub) — use the quay.io mirror.
FASTQC_IMAGE = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

FASTQC_THREADS = 2
DEFAULT_RESOURCES = flyte.Resources(cpu=FASTQC_THREADS, memory="6Gi")

fastqc = shell.create(
    name="fastqc",
    image=FASTQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        # FastQC allocates --memory per thread: two thirds of each CPU's share, at most 10000.
        FASTQC_MEMORY_MB=$(( TASK_MEMORY_MB * 2 / 3 / TASK_CPUS ))
        if [ $FASTQC_MEMORY_MB -gt 10000 ]; then FASTQC_MEMORY_MB=10000; fi
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        if [ ${{#R2[@]}} -eq 0 ]; then
            ln -s {{inputs.reads_1}} "$P.gz"
            FILES=("$P.gz")
        else
            ln -s {{inputs.reads_1}} "${{P}}_1.gz"
            ln -s "${{R2[0]}}" "${{P}}_2.gz"
            FILES=("${{P}}_1.gz" "${{P}}_2.gz")
        fi
        fastqc $ARGS --threads $TASK_CPUS --memory $FASTQC_MEMORY_MB -o {{outputs.results}} "${{FILES[@]}}"
    """,
)


# As upstream's module with a list of reads: one file is linked as
# <prefix>.gz, several as <prefix>_1.gz, <prefix>_2.gz, ... List inputs are
# staged in file-name order, which for Illumina names (..._L001_R1, _L001_R2,
# _L002_R1, ...) matches the samplesheet's order.
fastqc_reads = shell.create(
    name="fastqc_reads",
    image=FASTQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        # FastQC allocates --memory per thread: two thirds of each CPU's share, at most 10000.
        FASTQC_MEMORY_MB=$(( TASK_MEMORY_MB * 2 / 3 / TASK_CPUS ))
        if [ $FASTQC_MEMORY_MB -gt 10000 ]; then FASTQC_MEMORY_MB=10000; fi
        READS=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        FILES=()
        if [ ${{#READS[@]}} -eq 1 ]; then
            ln -s "${{READS[0]}}" "$P.gz"; FILES=("$P.gz")
        else
            i=1
            for r in "${{READS[@]}}"; do ln -s "$r" "${{P}}_$i.gz"; FILES+=("${{P}}_$i.gz"); i=$((i + 1)); done
        fi
        fastqc $ARGS --threads $TASK_CPUS --memory $FASTQC_MEMORY_MB -o {{outputs.results}} "${{FILES[@]}}"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "fastqc",
    fastqc.as_task(),
    fastqc_reads.as_task(),
)
