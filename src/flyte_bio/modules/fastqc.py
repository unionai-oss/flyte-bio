"""fastqc — read quality-control reports.

Exposes :data:`fastqc`, which writes FastQC's HTML report and data zip for
one sample's reads (one or two mates) into an output Dir.

As upstream, the reads are first symlinked to ``<prefix>.gz`` (single-end)
or ``<prefix>_1.gz`` / ``<prefix>_2.gz`` (paired) so the reports are named
after the sample, e.g. ``<prefix>_1_fastqc.html``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/fastqc:...` (Docker Hub) — use the quay.io mirror.
FASTQC_IMAGE = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

# FastQC allocates --memory per thread; keep threads x memory within the pod.
FASTQC_THREADS = 2
FASTQC_MEMORY_MB = 2048
DEFAULT_RESOURCES = flyte.Resources(cpu=FASTQC_THREADS, memory="6Gi")

fastqc = shell.create(
    name="fastqc",
    image=FASTQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
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
        fastqc $ARGS --threads {FASTQC_THREADS} --memory {FASTQC_MEMORY_MB} -o {{outputs.results}} "${{FILES[@]}}"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "fastqc",
    fastqc.as_task(),
)
