"""stardist — nucleus segmentation of 2D images with StarDist.

Exposes :func:`stardist`, which segments a 2D image (TIFF) with StarDist's
``stardist-predict2d`` and returns the label mask
``<image name without extension>.stardist.tif``. The model is a model folder
(``model_dir``) or, without one, a named pretrained model (``model``, e.g.
``2D_versatile_fluo``), which TensorFlow downloads at run time into a Keras
cache in the task's working directory, as ``HOME`` in the container may not
be writable. StarDist's options (``--n_tiles``, ``--prob_thresh``, ...) go in
``args``.

Upstream labels it process_gpu; the image carries CUDA libraries, and
TensorFlow uses a GPU when the task has one (request one with
``stardist_cmd.as_task().override(resources=...)``). On CPU, TensorFlow's
threads follow the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

STARDIST_IMAGE = "community.wave.seqera.io/library/python_pip_imagecodecs_nvidia-cublas-cu12_pruned:b668bcb6d531d350"
# Upstream labels it process_medium; large images need more memory (or --n_tiles in args).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

stardist_cmd = shell.create(
    name="stardist",
    image=STARDIST_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"image": File, "model": str, "model_dir": Dir | None, "args": str},
    defaults={"model": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGE=({{inputs.image}})
        MODEL={{inputs.model}}
        MODEL_DIR={{inputs.model_dir}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib" KERAS_HOME="$W/.keras"
        export OMP_NUM_THREADS=$TASK_CPUS TF_NUM_INTRAOP_THREADS=$TASK_CPUS TF_NUM_INTEROP_THREADS=$TASK_CPUS
        mkdir -p "$MPLCONFIGDIR" "$KERAS_HOME"
        EXTRA=()
        if [ -d "$MODEL_DIR" ]; then EXTRA+=(-m "$MODEL_DIR"); elif [ -n "$MODEL" ]; then EXTRA+=(-m "$MODEL"); fi
        stardist-predict2d \
            -i "${{IMAGE[0]}}" \
            -o "$OUT" \
            "${{EXTRA[@]}}" \
            $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("stardist", stardist_cmd.as_task())


async def stardist(image: File, model: str = "", model_dir: Dir | None = None, args: str = "") -> File:
    """The StarDist label mask of ``image``: ``<image name without extension>.stardist.tif``."""
    results = await stardist_cmd(image=image, model=model, model_dir=model_dir, args=args)
    files = await results.list_files()
    assert len(files) == 1, f"stardist wrote {len(files)} masks, expected 1"
    return files[0]
