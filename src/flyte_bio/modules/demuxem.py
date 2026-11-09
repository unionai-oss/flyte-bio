"""demuxem — demultiplex cell-hashing data with demuxEM.

Exposes :func:`demuxem`, which reads a 10x raw RNA count matrix (``.h5``)
and the matching hashtag count table (CSV, one row per hashtag, one column
per cell barcode) and assigns each droplet to a sample, a doublet or nothing
with demuxEM's expectation-maximisation. Outputs keep upstream's names:

- ``<prefix>_demux.zarr.zip``: the RNA matrix with the assignments;
- ``<prefix>.out.demuxEM.zarr.zip``: the hashtag matrix with demuxEM's
  statistics;
- with ``generate_diagnostic_plots`` or ``generate_gender_plot``, demuxEM's
  PDFs, kept in the result's ``results`` Dir.

``args`` go to demuxEM (e.g. ``--threads 1 --random-state 0``); ``genome``
is its ``--genome``. Upstream passes ``--generate-diagnostic-plots`` a value,
which demuxEM's option doesn't take; here it's a plain flag.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

DEMUXEM_IMAGE = "quay.io/biocontainers/demuxem:0.1.7.post1--pyhdfd78af_0"
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

demuxem_cmd = shell.create(
    name="demuxem",
    image=DEMUXEM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "rna_h5": File,
        "hto_csv": File,
        "prefix": str,
        "genome": str,
        "generate_gender_plot": str,
        "generate_diagnostic_plots": bool,
        "args": str,
    },
    defaults={"genome": "", "generate_gender_plot": "", "generate_diagnostic_plots": False, "args": ""},
    outputs={"results": Dir},
    script=r"""
        H5=({inputs.rna_h5}); CSV=({inputs.hto_csv})
        P={inputs.prefix}
        GENOME={inputs.genome}
        GENDER={inputs.generate_gender_plot}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        OPTS=()
        if [ -n "$GENDER" ]; then OPTS+=(--generate-gender-plot "$GENDER"); fi
        if [ -n "$GENOME" ]; then OPTS+=(--genome "$GENOME"); fi
        if [ {inputs.generate_diagnostic_plots} = true ]; then OPTS+=(--generate-diagnostic-plots); fi
        demuxEM $ARGS "${OPTS[@]}" "${H5[0]}" "${CSV[0]}" "$P" || exit 1
        mv "$P"* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("demuxem", demuxem_cmd.as_task())


@dataclass
class DemuxemResult:
    zarr: File  # <prefix>_demux.zarr.zip
    out_zarr: File  # <prefix>.out.demuxEM.zarr.zip
    results: Dir  # everything demuxEM wrote, including any plots


async def demuxem(
    rna_h5: File,
    hto_csv: File,
    prefix: str,
    genome: str = "",
    generate_gender_plot: str = "",
    generate_diagnostic_plots: bool = False,
    args: str = "",
) -> DemuxemResult:
    """demuxEM assignments for a raw 10x RNA matrix (``.h5``) and its hashtag count table (CSV)."""
    results = await demuxem_cmd(
        rna_h5=rna_h5,
        hto_csv=hto_csv,
        prefix=prefix,
        genome=genome,
        generate_gender_plot=generate_gender_plot,
        generate_diagnostic_plots=generate_diagnostic_plots,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"demuxem wrote no {name}")
        return found

    return DemuxemResult(
        zarr=await pick(f"{prefix}_demux.zarr.zip"),
        out_zarr=await pick(f"{prefix}.out.demuxEM.zarr.zip"),
        results=results,
    )
