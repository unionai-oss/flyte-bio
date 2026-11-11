"""numorph — light-sheet microscopy preprocessing and nuclei counting with NuMorph.

Exposes three functions, one per NuMorph stage. Each takes its image stack as
an archive (``.zip``, ``.tar.gz``, ``.tgz`` or ``.tar``) of the image
directory rather than a Dir, as stacks hold hundreds of TIFFs, more than a
task's input staging copes with. An archive holding a single top-level
directory is read from inside that directory.

- :func:`numorph_intensity` runs NuMorph's ``intensity`` preprocessing stage
  (shading correction and intensity normalisation between tile stacks) on raw
  tiles, configured by a NuMorph parameter CSV. It returns, under upstream's
  names, ``variables/`` (``adj_params.mat``, ``path_table.mat``,
  ``thresholds.mat``), ``samples/`` (``intensity_adjustment/*.{tif,png}``)
  and ``NM_variables.mat``.
- :func:`numorph_resample` runs NuMorph's ``resample`` analysis stage on
  stitched images (one tile per channel), returning ``resampled/``:
  ``<sample>_C<n>_<marker>_<resolution>.nii`` per channel.
- :func:`numorph_3dunet` counts cell nuclei in stitched images with NuMorph's
  3D U-Net and a trained model (``.h5``), returning ``<sample>.csv`` (the
  centroids) and ``<sample>_counts.csv``. It uses the GPU when the task has
  one (request one with ``numorph_3dunet_cmd.override(...)``) and
  otherwise runs on CPU. ``args`` passes extra options (``--chunk_size``,
  ``--pred_threshold``, ...). As upstream, it takes no parameter CSV.

Licences: NuMorph and its images are MIT. ``intensity`` and ``resample`` are
compiled MATLAB applications, so their images ship the MATLAB Runtime R2023a,
which MathWorks licenses under its MATLAB Runtime licence (free to use and
redistribute, not open source). MATLAB chooses its own thread count. The 3D
U-Net image is built on NVIDIA's TensorFlow 1.15 container (20.10-tf1), which
carries NVIDIA's Deep Learning Container licence; its TensorFlow needs a CPU
with AVX.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PREPROCESSING_IMAGE = "quay.io/nf-core/numorph_preprocessing:1.0.0"
ANALYZE_IMAGE = "quay.io/nf-core/numorph_analyze:1.0.0"
UNET_IMAGE = "quay.io/nf-core/numorph-3dunet:1.0.9"
# Upstream labels intensity and resample process_low. Whole-brain stacks are
# tens to hundreds of GB and need far more memory and disk than these defaults.
MATLAB_RESOURCES = flyte.Resources(cpu=2, memory="8Gi", disk="20Gi")
# Upstream labels it process_gpu.
UNET_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")

# Extracts the archive IMAGES[0] into $IMAGE_DIR, descending into a lone top-level directory.
EXTRACT_IMAGES = r"""
        IMAGE_DIR=$(mktemp -d)
        case "${IMAGES[0]}" in
            *.zip) unzip -q "${IMAGES[0]}" -d "$IMAGE_DIR" || exit 1 ;;
            *) tar -xf "${IMAGES[0]}" -C "$IMAGE_DIR" || exit 1 ;;
        esac
        ENTRIES=("$IMAGE_DIR"/*)
        if [ ${#ENTRIES[@]} -eq 1 ] && [ -d "${ENTRIES[0]}" ]; then IMAGE_DIR=${ENTRIES[0]}; fi
""".strip()

# The MATLAB Runtime unpacks the application into its cache.
MATLAB_SETUP = r"""
        W=$(mktemp -d); cd "$W"
        export MCR_CACHE_ROOT="$W/.mcr" HOME="$W"
        mkdir -p "$MCR_CACHE_ROOT"
""".strip()

numorph_intensity_cmd = shell.create(
    name="numorph_intensity",
    image=PREPROCESSING_IMAGE,
    resources=MATLAB_RESOURCES,
    inputs={"images": File, "parameter_file": File, "sample_name": str},
    outputs={"results": Dir},
    script=rf"""
        IMAGES=({{inputs.images}})
        PARAMS=({{inputs.parameter_file}})
        SAMPLE={{inputs.sample_name}}
        OUT={{outputs.results}}
        {MATLAB_SETUP}
        {EXTRACT_IMAGES}
        mkdir -p results
        numorph_preprocessing 'input_dir' "$IMAGE_DIR" 'output_dir' "$W/results" \
            'parameter_file' "${{PARAMS[0]}}" 'sample_name' "$SAMPLE" 'stage' 'intensity' || exit 1
        for P in variables samples NM_variables.mat; do
            if [ ! -e "results/$P" ]; then echo "numorph_preprocessing wrote no results/$P" >&2; exit 1; fi
            mv "results/$P" "$OUT"/
        done
    """,
)

numorph_resample_cmd = shell.create(
    name="numorph_resample",
    image=ANALYZE_IMAGE,
    resources=MATLAB_RESOURCES,
    inputs={"images": File, "parameter_file": File, "sample_name": str},
    outputs={"results": Dir},
    script=rf"""
        IMAGES=({{inputs.images}})
        PARAMS=({{inputs.parameter_file}})
        SAMPLE={{inputs.sample_name}}
        OUT={{outputs.results}}
        {MATLAB_SETUP}
        {EXTRACT_IMAGES}
        mkdir -p results/stitched
        ln -s "$IMAGE_DIR"/* results/stitched/
        numorph_analyze 'input_dir' "$W/results/stitched" \
            'output_dir' "$W/results" \
            'parameter_file' "${{PARAMS[0]}}" \
            'sample_name' "$SAMPLE" \
            'stage' 'resample' \
            'NM_variables' '' \
            'use_processed_images' 'stitched' || exit 1
        shopt -s nullglob
        RESAMPLED=(results/resampled/*)
        if [ ${{#RESAMPLED[@]}} -eq 0 ]; then echo "numorph_analyze wrote nothing to results/resampled" >&2; exit 1; fi
        mv "${{RESAMPLED[@]}}" "$OUT"/
    """,
)

numorph_3dunet_cmd = shell.create(
    name="numorph_3dunet",
    image=UNET_IMAGE,
    resources=UNET_RESOURCES,
    inputs={"images": File, "model_file": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IMAGES=({{inputs.images}})
        MODEL=({{inputs.model_file}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" OMP_NUM_THREADS=$TASK_CPUS MKL_NUM_THREADS=$TASK_CPUS
        if command -v nvidia-smi > /dev/null && nvidia-smi -L > /dev/null 2>&1; then
            echo "Checking GPU access:"
            nvidia-smi
        fi
        {EXTRACT_IMAGES}
        mkdir -p "$P"
        numorph_3dunet \
            -i "$IMAGE_DIR" \
            -o "$P" \
            --model_file "${{MODEL[0]}}" \
            --sample_id "$P" \
            $ARGS || exit 1
        shopt -s nullglob
        COUNTS=("$P"/*)
        if [ ${{#COUNTS[@]}} -eq 0 ]; then echo "numorph_3dunet wrote nothing" >&2; exit 1; fi
        mv "${{COUNTS[@]}}" "$OUT"/
    """,
)

numorph_preprocessing_env = flyte.TaskEnvironment.from_task("numorph_preprocessing", numorph_intensity_cmd.as_task())
numorph_analyze_env = flyte.TaskEnvironment.from_task("numorph_analyze", numorph_resample_cmd.as_task())
numorph_3dunet_env = flyte.TaskEnvironment.from_task("numorph_3dunet", numorph_3dunet_cmd.as_task())
env = flyte.TaskEnvironment(
    name="numorph",
    depends_on=[numorph_preprocessing_env, numorph_analyze_env, numorph_3dunet_env],
)


@dataclass
class NumorphIntensityResult:
    variables: list[File]  # variables/*.mat, sorted by name
    samples: list[File]  # samples/**, sorted by path under samples/
    nm_variables: File  # NM_variables.mat
    results: Dir  # all of the above


@dataclass
class NumorphUnetResult:
    centroids: File  # <prefix>.csv
    counts: File  # <prefix>_counts.csv
    results: Dir  # all files numorph_3dunet wrote


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def numorph_intensity(images: File, parameter_file: File, sample_name: str) -> NumorphIntensityResult:
    """Intensity-adjust the raw tiles in the archive ``images`` as ``parameter_file`` sets."""
    results = await numorph_intensity_cmd(images=images, parameter_file=parameter_file, sample_name=sample_name)
    base = results.path.rstrip("/")
    paths = sorted([f.path[len(base) :].lstrip("/") async for f in results.walk()])

    def under(prefix: str) -> list[File]:
        return [File.from_existing_remote(f"{base}/{p}") for p in paths if p.startswith(prefix)]

    return NumorphIntensityResult(
        variables=under("variables/"),
        samples=under("samples/"),
        nm_variables=File.from_existing_remote(f"{base}/NM_variables.mat"),
        results=results,
    )


async def numorph_resample(images: File, parameter_file: File, sample_name: str) -> list[File]:
    """Resample the stitched images in the archive ``images``: ``resampled/*``, sorted by name."""
    results = await numorph_resample_cmd(images=images, parameter_file=parameter_file, sample_name=sample_name)
    return sorted(await results.list_files(), key=file_name)


async def numorph_3dunet(images: File, model_file: File, prefix: str, args: str = "") -> NumorphUnetResult:
    """Count the nuclei in the stitched images in the archive ``images`` with the 3D U-Net ``model_file``."""
    results = await numorph_3dunet_cmd(images=images, model_file=model_file, prefix=prefix, args=args)
    names = {file_name(f): f for f in await results.list_files()}
    for wanted in (f"{prefix}.csv", f"{prefix}_counts.csv"):
        if wanted not in names:
            raise FileNotFoundError(f"numorph_3dunet wrote no {wanted} (wrote {sorted(names)})")
    return NumorphUnetResult(centroids=names[f"{prefix}.csv"], counts=names[f"{prefix}_counts.csv"], results=results)
