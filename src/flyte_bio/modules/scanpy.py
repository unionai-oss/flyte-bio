"""scanpy — AnnData conversions of single-cell count matrices.

Exposes:

- :func:`mtx_to_h5ad_simpleaf` — a sample's simpleaf quantification as
  ``<sample>_<input_type>_matrix.h5ad`` (cells and genes sorted, gene IDs
  without version suffixes, the sample recorded in ``obs``).
- :func:`concat_h5ad` — those per-sample files concatenated into
  ``<prefix>_matrix.h5ad``, with each cell annotated from the samplesheet.

Both run upstream's scripts (vendored) on upstream's scanpy image.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

SCANPY_IMAGE = "community.wave.seqera.io/library/scanpy:1.10.2--e83da2205b92a538"

DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# The scripts point numba's cache at the working directory, so run in a scratch one.
mtx_to_h5ad_simpleaf_cmd = shell.create(
    name="mtx_to_h5ad_simpleaf",
    image=SCANPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "quant": Dir, "sample": str, "input_type": str},
    outputs={"results": Dir},
    script=r"""
        S=({inputs.script})
        QUANT={inputs.quant}
        W=$(mktemp -d); cd "$W"
        python "${S[0]}" --input "$QUANT" --sample {inputs.sample} \
            --output {outputs.results}/{inputs.sample}_{inputs.input_type}_matrix.h5ad
    """,
)

# Upstream stages the h5ads flat in the working directory and keys each by its
# path there, so link them into a scratch directory and run from it.
concat_h5ad_cmd = shell.create(
    name="concat_h5ad",
    image=SCANPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "h5ads": list[File], "samplesheet": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        S=({inputs.script}); H=({inputs.h5ads}); SHEET=({inputs.samplesheet})
        W=$(mktemp -d); cd "$W"
        for f in "${H[@]}"; do ln -s "$f" "$(basename "$f")"; done
        python "${S[0]}" --samplesheet "${SHEET[0]}" --output {outputs.results}/{inputs.prefix}_matrix.h5ad
    """,
)


env = flyte.TaskEnvironment.from_task(
    "scanpy",
    mtx_to_h5ad_simpleaf_cmd.as_task(),
    concat_h5ad_cmd.as_task(),
)


async def pick(results: Dir, name: str, task: str) -> File:
    f = await results.get_file(name)
    if f is None:
        raise FileNotFoundError(f"{task} wrote no {name}")
    return f


async def mtx_to_h5ad_simpleaf(quant: Dir, sample: str, input_type: str = "raw") -> File:
    """``<sample>_<input_type>_matrix.h5ad`` from a simpleaf ``af_quant`` directory."""
    script = await File.from_local(str(path("mtx_to_h5ad_simpleaf.py")))
    results = await mtx_to_h5ad_simpleaf_cmd(script=script, quant=quant, sample=sample, input_type=input_type)
    return await pick(results, f"{sample}_{input_type}_matrix.h5ad", "mtx_to_h5ad_simpleaf")


async def concat_h5ad(h5ads: list[File], samplesheet: File, prefix: str = "combined_raw") -> File:
    """``<prefix>_matrix.h5ad``: the per-sample ``*_matrix.h5ad`` files concatenated, annotated from the samplesheet.

    Each cell's ``sample`` is its file's name without ``_matrix.h5ad``, as upstream.
    """
    names = [h.path.rstrip("/").rsplit("/", 1)[-1] for h in h5ads]
    if len(set(names)) != len(names):
        raise ValueError(f"concat_h5ad needs distinct file names, got {names}")
    script = await File.from_local(str(path("concat_h5ad.py")))
    results = await concat_h5ad_cmd(script=script, h5ads=h5ads, samplesheet=samplesheet, prefix=prefix)
    return await pick(results, f"{prefix}_matrix.h5ad", "concat_h5ad")
