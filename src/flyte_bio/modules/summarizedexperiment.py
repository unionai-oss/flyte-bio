"""summarizedexperiment — bundle count matrices into a SummarizedExperiment.

Reads one or more feature-by-sample matrices as assays, optionally attaches
row (feature) and column (sample) metadata, and saves the object as RDS. The
R helper is vendored and staged into the stock SummarizedExperiment
biocontainer as a File input.

Assay names are passed explicitly as ``name=filename`` pairs rather than by
position: a ``list[File]`` is staged under original basenames and globbed
alphabetically, so file order can't carry the names. The helper picks each
file's delimiter from its extension, so matrices and metadata must be
``.tsv``/``.txt``/``.csv`` files.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

SUMMARIZEDEXPERIMENT_IMAGE = "quay.io/biocontainers/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

summarized_experiment_cmd = shell.create(
    name="summarized_experiment",
    image=SUMMARIZEDEXPERIMENT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "script": File,
        "matrices": list[File],
        "assays": str,
        "rowdata": File | None,
        "coldata": File | None,
        "prefix": str,
    },
    outputs={"results": Dir},
    script=r"""
        M=({inputs.matrices})
        shopt -s nullglob; ROW=({inputs.rowdata}); COL=({inputs.coldata}); shopt -u nullglob
        OPTS=()
        [ ${#ROW[@]} -gt 0 ] && OPTS+=(--rowdata "${ROW[0]}")
        [ ${#COL[@]} -gt 0 ] && OPTS+=(--coldata "${COL[0]}")
        Rscript {inputs.script} \
            --matrices "$(dirname "${M[0]}")" \
            --assays {inputs.assays} \
            --prefix {outputs.results}/{inputs.prefix} \
            "${OPTS[@]}"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "summarizedexperiment",
    summarized_experiment_cmd.as_task(),
)


async def summarized_experiment(
    assays: dict[str, File],
    rowdata: File | None = None,
    coldata: File | None = None,
    prefix: str = "all_samples",
) -> File:
    """Return the ``<prefix>.SummarizedExperiment.rds`` built from ``assays`` (name → matrix)."""
    if not assays:
        raise ValueError("summarized_experiment needs at least one assay")
    filenames = {name: f.path.rstrip("/").rsplit("/", 1)[-1] for name, f in assays.items()}
    if len(set(filenames.values())) != len(filenames):
        raise ValueError(f"assay files must have distinct basenames: {filenames}")
    for name in filenames:
        if "=" in name or "," in name:
            raise ValueError(f"assay name {name!r} must not contain '=' or ','")

    script = await File.from_local(str(path("summarizedexperiment.r")))
    results = await summarized_experiment_cmd(
        script=script,
        matrices=list(assays.values()),
        assays=",".join(f"{name}={fname}" for name, fname in filenames.items()),
        rowdata=rowdata,
        coldata=coldata,
        prefix=prefix,
    )
    rds = await results.get_file(f"{prefix}.SummarizedExperiment.rds")
    if rds is None:
        raise FileNotFoundError(f"summarized_experiment wrote no {prefix}.SummarizedExperiment.rds")
    return rds
