"""cellpose — cell segmentation of microscopy images with Cellpose.

Exposes :func:`cellpose`, which segments one image (TIFF, OME-TIFF, PNG, ...)
with Cellpose's command line (``--save_tif --verbose``, as upstream) and
returns what it wrote, under Cellpose's names:

- ``*masks.tif``: the label masks (``<image>_cp_masks.tif``, and with
  ``--save_flows`` in ``args`` also ``<image>_flows_cp_masks.tif`` and
  ``<image>_dP_cp_masks.tif``);
- ``*flows.tif``: flows, if Cellpose writes files so named;
- ``*seg.npy``: Cellpose's ``_seg.npy`` (unless ``args`` has ``--no_npy``).

Without ``model``, Cellpose uses its default pretrained model, downloading it
at run time (into the task's working directory, as ``HOME`` in the container
isn't writable). Pass a model file to use it with ``--pretrained_model``.
Cellpose uses the GPU (``--use_gpu``) when the task has one, as upstream when
the process has an accelerator; request one with
``cellpose_cmd.as_task().override(resources=...)``. On CPU, PyTorch's threads
follow the task's CPUs.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CELLPOSE_IMAGE = "community.wave.seqera.io/library/python_pip_cellpose:fdf7a8c3a305a26e"
# Upstream labels it process_medium and process_gpu; large images need more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

cellpose_cmd = shell.create(
    name="cellpose",
    image=CELLPOSE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "model": File | None, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGE=({{inputs.image}})
        shopt -s nullglob; MODEL=({{inputs.model}}); shopt -u nullglob
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export OMP_NUM_THREADS=$TASK_CPUS MKL_NUM_THREADS=$TASK_CPUS
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib" CELLPOSE_LOCAL_MODELS_PATH="$W/.cellpose"
        mkdir -p "$MPLCONFIGDIR" "$CELLPOSE_LOCAL_MODELS_PATH" images
        # Cellpose writes its outputs next to the image.
        ln -s "${{IMAGE[0]}}" images/
        EXTRA=()
        if command -v nvidia-smi > /dev/null && nvidia-smi -L > /dev/null 2>&1; then EXTRA+=(--use_gpu); fi
        if [ ${{#MODEL[@]}} -gt 0 ]; then EXTRA+=(--pretrained_model "${{MODEL[0]}}"); fi
        cellpose \
            --image_path "images/$(basename "${{IMAGE[0]}}")" \
            --save_tif \
            --verbose \
            "${{EXTRA[@]}}" \
            $ARGS
        shopt -s nullglob
        MASKS=(images/*masks.tif)
        if [ ${{#MASKS[@]}} -eq 0 ]; then echo "cellpose wrote no masks" >&2; exit 1; fi
        mv images/*masks.tif images/*flows.tif images/*seg.npy "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("cellpose", cellpose_cmd.as_task())


@dataclass
class CellposeResult:
    masks: list[File]  # *masks.tif, sorted by name
    flows: list[File]  # *flows.tif
    cells: list[File]  # *seg.npy
    results: Dir  # all of the above


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def cellpose(image: File, model: File | None = None, args: str = "") -> CellposeResult:
    """Segment the cells in ``image`` with Cellpose, using ``model`` or the default pretrained model."""
    results = await cellpose_cmd(image=image, model=model, args=args)
    files = sorted(await results.list_files(), key=file_name)

    def with_suffix(suffix: str) -> list[File]:
        return [f for f in files if file_name(f).endswith(suffix)]

    return CellposeResult(
        masks=with_suffix("masks.tif"), flows=with_suffix("flows.tif"), cells=with_suffix("seg.npy"), results=results
    )
