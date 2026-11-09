"""learnmsa — protein multiple sequence alignment with learnMSA.

Exposes :func:`learnmsa_align`, which aligns the sequences in ``fasta`` with
learnMSA's profile hidden Markov models and writes ``<prefix>.aln``.
learnMSA's options go in ``args``.

learnMSA trains with TensorFlow, which uses a GPU when it finds one. These
tasks run on CPUs only: no GPU is requested, and ``CUDA_VISIBLE_DEVICES`` is
emptied so TensorFlow doesn't look for one. Training is randomised, so the
alignment can differ between runs (upstream's test treats it as unstable).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

LEARNMSA_IMAGE = "quay.io/biocontainers/learnmsa:2.0.16--pyhdfd78af_0"
# Upstream labels this process_medium. Large families train faster on a GPU.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

learnmsa_align_cmd = shell.create(
    name="learnmsa_align",
    image=LEARNMSA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # CPU only; TensorFlow and matplotlib write caches under HOME.
        export CUDA_VISIBLE_DEVICES="" HOME="$W" MPLCONFIGDIR="$W/matplotlib"
        export TF_NUM_INTRAOP_THREADS=$TASK_CPUS TF_NUM_INTEROP_THREADS=$TASK_CPUS
        learnMSA \
            -i "${{FA[0]}}" \
            -o "$OUT/$P.aln" \
            $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("learnmsa", learnmsa_align_cmd.as_task())


async def learnmsa_align(fasta: File, prefix: str = "sample", args: str = "") -> File:
    """Align ``fasta`` with learnMSA: ``<prefix>.aln``."""
    results = await learnmsa_align_cmd(fasta=fasta, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.aln")
    if found is None:
        raise FileNotFoundError(f"learnmsa wrote no {prefix}.aln")
    return found
