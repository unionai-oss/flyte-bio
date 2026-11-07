"""deepcell — whole-cell and nuclear segmentation of multiplexed images with DeepCell's Mesmer.

Exposes :func:`deepcell_mesmer`, which segments a nuclear image (optionally
with a separate membrane image) with deepcell-applications' ``run_app.py
mesmer --squeeze`` and returns the label mask ``<prefix>.tif``. Mesmer's
options (``--image-mpp``, ``--compartment``, ``--nuclear-channel``,
``--membrane-channel``, ...) go in ``args``.

The image ships Mesmer's pretrained model under ``/root/.keras``. As ``HOME``
in the container may not be writable, the task links that copy into a
Keras cache in its working directory (where TensorFlow extracts it), or, if
it can't read it, lets TensorFlow download the model there. TensorFlow's
threads follow the task's CPUs. The image is CPU-only.

DeepCell and its models are licensed for non-commercial use only (a
modified Apache 2.0 licence).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MESMER_IMAGE = "quay.io/nf-core/deepcell_mesmer:0.4.1_noentry"
# Upstream labels it process_high; this is enough for test-sized images. Whole slides need more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

deepcell_mesmer_cmd = shell.create(
    name="deepcell_mesmer",
    image=MESMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "membrane_image": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGE=({{inputs.image}})
        shopt -s nullglob; MEMBRANE=({{inputs.membrane_image}}); shopt -u nullglob
        PREFIX={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib"
        export OMP_NUM_THREADS=$TASK_CPUS TF_NUM_INTRAOP_THREADS=$TASK_CPUS TF_NUM_INTEROP_THREADS=$TASK_CPUS
        mkdir -p "$MPLCONFIGDIR" "$W/.keras/models"
        BUNDLED=/root/.keras/models/MultiplexSegmentation.tgz
        if [ -r "$BUNDLED" ]; then ln -s "$BUNDLED" "$W/.keras/models/"; fi
        EXTRA=()
        if [ ${{#MEMBRANE[@]}} -gt 0 ]; then EXTRA+=(--membrane-image "${{MEMBRANE[0]}}"); fi
        python /usr/src/app/run_app.py mesmer \
            --squeeze \
            --nuclear-image "${{IMAGE[0]}}" \
            --output-directory "$OUT" \
            --output-name "$PREFIX.tif" \
            "${{EXTRA[@]}}" \
            $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("deepcell", deepcell_mesmer_cmd.as_task())


async def deepcell_mesmer(image: File, prefix: str, membrane_image: File | None = None, args: str = "") -> File:
    """The Mesmer label mask of ``image`` (and ``membrane_image``, if given): ``<prefix>.tif``."""
    results = await deepcell_mesmer_cmd(image=image, membrane_image=membrane_image, prefix=prefix, args=args)
    return await results.get_file(f"{prefix}.tif")
