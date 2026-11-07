"""s4pred — protein secondary structure prediction with S4PRED.

Exposes :func:`s4pred_runmodel`, which predicts the secondary structure of
every sequence in a protein FASTA with S4PRED's ``run_model``. As upstream it
saves one PSIPRED-style ``.ss2`` file per sequence (``--save-files``), named
after the sequence's ID, and returns them in a Dir (upstream's
``<prefix>`` folder). ``run_model`` runs with ``--threads`` set to the task's
CPUs. The model weights ship in the container.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

S4PRED_IMAGE = "quay.io/biocontainers/s4pred:1.2.1--pyhdfd78af_1"
# Upstream: process_high (12 CPUs, 72 GB). The test proteome runs in far less;
# raise both for large proteomes.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

s4pred_runmodel_cmd = shell.create(
    name="s4pred_runmodel",
    image=S4PRED_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"preds": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        ARGS={{inputs.args}}
        OUT={{outputs.preds}}
        W=$(mktemp -d); cd "$W"

        run_model \
            $ARGS \
            --threads $TASK_CPUS \
            --save-files \
            --outdir "$OUT" \
            "${{FA[0]}}" > /dev/null || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("s4pred", s4pred_runmodel_cmd.as_task())


async def s4pred_runmodel(fasta: File, args: str = "") -> Dir:
    """Secondary structure of each protein in ``fasta``: a Dir of ``<sequence id>.ss2``."""
    return await s4pred_runmodel_cmd(fasta=fasta, args=args)
