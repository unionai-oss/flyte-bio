"""cellbender — remove ambient RNA and empty droplets from a raw count matrix.

Exposes :func:`cellbender_removebackground`, which runs ``cellbender
remove-background`` on an unfiltered h5ad and writes upstream's outputs:
``<prefix>.h5`` (the denoised counts), ``<prefix>_filtered.h5``,
``<prefix>_posterior.h5``, ``<prefix>_cell_barcodes.csv``,
``<prefix>_metrics.csv``, ``<prefix>_report.html`` (when CellBender writes
one), ``<prefix>.pdf``, ``<prefix>.log`` and the training checkpoint
``ckpt.tar.gz``.

Like upstream, it runs on the CPU image by default, and with ``gpu=True``
on Broad's CUDA image with ``--cuda`` and one GPU.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CELLBENDER_IMAGE = "community.wave.seqera.io/library/cellbender_webcolors:156d413fdfc16cdb"
CELLBENDER_GPU_IMAGE = "us.gcr.io/broad-dsde-methods/cellbender:0.3.2"

DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
GPU_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", gpu=1)

# CellBender writes its checkpoint to the working directory and its temporary
# files under TMPDIR, so run in a scratch directory and move the outputs.
SCRIPT = f"""
    {TASK_LIMITS}
""" + r"""
    H=({inputs.h5ad})
    ARGS={inputs.args}
    W=$(mktemp -d); cd "$W"
    TMPDIR=. cellbender remove-background $ARGS --cpu-threads $TASK_CPUS --estimator-multiple-cpu \
        CUDA_FLAG --input "${H[0]}" --output {inputs.prefix}.h5
    mv {inputs.prefix}.h5 {inputs.prefix}_filtered.h5 {inputs.prefix}_posterior.h5 \
        {inputs.prefix}_cell_barcodes.csv {inputs.prefix}_metrics.csv {inputs.prefix}.pdf \
        {inputs.prefix}.log ckpt.tar.gz {outputs.results}/
    if [ -f {inputs.prefix}_report.html ]; then mv {inputs.prefix}_report.html {outputs.results}/; fi
"""

INPUTS = {"h5ad": File, "prefix": str, "args": str}

cellbender_removebackground_cmd = shell.create(
    name="cellbender_removebackground",
    image=CELLBENDER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs=INPUTS,
    defaults={"args": ""},
    outputs={"results": Dir},
    script=SCRIPT.replace("CUDA_FLAG", ""),
)

cellbender_removebackground_gpu_cmd = shell.create(
    name="cellbender_removebackground_gpu",
    image=CELLBENDER_GPU_IMAGE,
    resources=GPU_RESOURCES,
    inputs=INPUTS,
    defaults={"args": ""},
    outputs={"results": Dir},
    script=SCRIPT.replace("CUDA_FLAG", "--cuda"),
)


# One image per task, so one env each; `env` aggregates them.
cpu_env = flyte.TaskEnvironment.from_task("cellbender_cpu", cellbender_removebackground_cmd.as_task())
gpu_env = flyte.TaskEnvironment.from_task("cellbender_gpu", cellbender_removebackground_gpu_cmd.as_task())
env = flyte.TaskEnvironment(name="cellbender", depends_on=[cpu_env, gpu_env])


@dataclass
class CellbenderResult:
    h5: File  # <prefix>.h5: denoised counts of every droplet
    filtered_h5: File  # <prefix>_filtered.h5: denoised counts of the called cells
    posterior_h5: File  # <prefix>_posterior.h5
    barcodes: File  # <prefix>_cell_barcodes.csv: the called cells, headerless
    metrics: File  # <prefix>_metrics.csv
    report: File | None  # <prefix>_report.html
    pdf: File  # <prefix>.pdf
    log: File  # <prefix>.log
    checkpoint: File  # ckpt.tar.gz


async def cellbender_removebackground(h5ad: File, prefix: str, gpu: bool = False, args: str = "") -> CellbenderResult:
    """``cellbender remove-background`` on an unfiltered h5ad; ``args`` go to CellBender (e.g. ``--epochs 150``)."""
    cmd = cellbender_removebackground_gpu_cmd if gpu else cellbender_removebackground_cmd
    results = await cmd(h5ad=h5ad, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        f = await results.get_file(name)
        if f is None:
            raise FileNotFoundError(f"cellbender_removebackground wrote no {name}")
        return f

    return CellbenderResult(
        h5=await pick(f"{prefix}.h5"),
        filtered_h5=await pick(f"{prefix}_filtered.h5"),
        posterior_h5=await pick(f"{prefix}_posterior.h5"),
        barcodes=await pick(f"{prefix}_cell_barcodes.csv"),
        metrics=await pick(f"{prefix}_metrics.csv"),
        report=await results.get_file(f"{prefix}_report.html"),
        pdf=await pick(f"{prefix}.pdf"),
        log=await pick(f"{prefix}.log"),
        checkpoint=await pick("ckpt.tar.gz"),
    )
