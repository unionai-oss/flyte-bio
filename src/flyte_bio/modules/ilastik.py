"""ilastik — interactive-learning image analysis, run headless with a trained project.

Exposes:

- :func:`ilastik_pixelclassification`: run a trained Pixel Classification
  project (``.ilp``) on an image, exporting ``export_source`` (``probabilities``,
  ``simple segmentation``, ``uncertainty``, ``features`` or ``labels``) as
  ``<prefix>.<output_format>``. ``output_format`` must match the export format
  saved in the project (e.g. ``h5``).
- :func:`ilastik_multicut`: run a trained Boundary-based Segmentation with
  Multicut project on an image and its boundary probabilities, exporting the
  multicut segmentation as ``<prefix>.tiff``.

Both run ilastik read-only on the project (``--headless --readonly 1``), as
upstream, and pass ``args`` on to ilastik. ilastik's threads
(``LAZYFLOW_THREADS``) follow the task's CPUs and its memory budget
(``LAZYFLOW_TOTAL_RAM_MB``) is 80% of the task's memory.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's image is Docker Hub's biocontainers/ilastik, which has no quay.io mirror (it isn't a bioconda build).
ILASTIK_IMAGE = "docker.io/biocontainers/ilastik:1.4.0_cv1"
RUN_ILASTIK = "/opt/ilastik-1.4.0-Linux/run_ilastik.sh"
# Upstream labels pixel classification process_single and multicut process_low; large images need more.
PIXELCLASSIFICATION_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
MULTICUT_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")

ILASTIK_SETUP = rf"""
        {TASK_LIMITS}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W" XDG_CACHE_HOME="$W/.cache"
        export LAZYFLOW_THREADS=$TASK_CPUS LAZYFLOW_TOTAL_RAM_MB=$(( TASK_MEMORY_MB * 8 / 10 ))
"""

ilastik_pixelclassification_cmd = shell.create(
    name="ilastik_pixelclassification",
    image=ILASTIK_IMAGE,
    resources=PIXELCLASSIFICATION_RESOURCES,
    inputs={
        "image": File,
        "project": File,
        "output_format": str,
        "export_source": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        IMAGE=({{inputs.image}})
        PROJECT=({{inputs.project}})
        FORMAT={{inputs.output_format}}
        SOURCE={{inputs.export_source}}
        PREFIX={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {ILASTIK_SETUP}
        {RUN_ILASTIK} \
            --headless \
            --readonly 1 \
            --project="${{PROJECT[0]}}" \
            --output_filename_format="$OUT/$PREFIX.$FORMAT" \
            --export_source="$SOURCE" \
            $ARGS \
            "${{IMAGE[0]}}" || exit 1
        [ -s "$OUT/$PREFIX.$FORMAT" ] || {{ echo "ilastik wrote no $PREFIX.$FORMAT" >&2; exit 1; }}
    """,
)

ilastik_multicut_cmd = shell.create(
    name="ilastik_multicut",
    image=ILASTIK_IMAGE,
    resources=MULTICUT_RESOURCES,
    inputs={"image": File, "probabilities": File, "project": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        IMAGE=({{inputs.image}})
        PROBABILITIES=({{inputs.probabilities}})
        PROJECT=({{inputs.project}})
        PREFIX={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {ILASTIK_SETUP}
        {RUN_ILASTIK} \
            --headless \
            --readonly 1 \
            --project="${{PROJECT[0]}}" \
            --raw_data="${{IMAGE[0]}}" \
            --probabilities="${{PROBABILITIES[0]}}" \
            --export_source="Multicut Segmentation" \
            --output_filename_format="$OUT/$PREFIX.tiff" \
            $ARGS || exit 1
        [ -s "$OUT/$PREFIX.tiff" ] || {{ echo "ilastik wrote no $PREFIX.tiff" >&2; exit 1; }}
    """,
)

env = flyte.TaskEnvironment.from_task(
    "ilastik", ilastik_pixelclassification_cmd.as_task(), ilastik_multicut_cmd.as_task()
)


async def ilastik_pixelclassification(
    image: File, project: File, output_format: str, export_source: str, prefix: str, args: str = ""
) -> File:
    """``export_source`` of the Pixel Classification ``project`` run on ``image``: ``<prefix>.<output_format>``."""
    results = await ilastik_pixelclassification_cmd(
        image=image,
        project=project,
        output_format=output_format,
        export_source=export_source,
        prefix=prefix,
        args=args,
    )
    return await results.get_file(f"{prefix}.{output_format}")


async def ilastik_multicut(image: File, probabilities: File, project: File, prefix: str, args: str = "") -> File:
    """The multicut segmentation of ``image`` and its boundary ``probabilities`` by ``project``: ``<prefix>.tiff``."""
    results = await ilastik_multicut_cmd(
        image=image, probabilities=probabilities, project=project, prefix=prefix, args=args
    )
    return await results.get_file(f"{prefix}.tiff")
