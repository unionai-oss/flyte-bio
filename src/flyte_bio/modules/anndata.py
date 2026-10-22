"""anndata — subset an AnnData h5ad.

Exposes :func:`anndata_barcodes`, which keeps an h5ad's cells whose
barcodes are listed in a CSV and writes ``<prefix>.h5ad``, running
upstream's script (vendored) on upstream's image.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

ANNDATA_IMAGE = "community.wave.seqera.io/library/anndata_pandas_python_pyyaml:3f88bccf1cdade40"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

anndata_barcodes_cmd = shell.create(
    name="anndata_barcodes",
    image=ANNDATA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "h5ad": File, "barcodes": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        S=({inputs.script}); H=({inputs.h5ad}); B=({inputs.barcodes})
        python3 "${S[0]}" --h5ad "${H[0]}" --barcodes "${B[0]}" --output {outputs.results}/{inputs.prefix}.h5ad
    """,
)


env = flyte.TaskEnvironment.from_task(
    "anndata",
    anndata_barcodes_cmd.as_task(),
)


async def anndata_barcodes(h5ad: File, barcodes: File, prefix: str) -> File:
    """``<prefix>.h5ad``: ``h5ad``'s cells listed in ``barcodes`` (headerless CSV), in its order."""
    script = await File.from_local(str(path("anndata_barcodes.py")))
    results = await anndata_barcodes_cmd(script=script, h5ad=h5ad, barcodes=barcodes, prefix=prefix)
    f = await results.get_file(f"{prefix}.h5ad")
    if f is None:
        raise FileNotFoundError(f"anndata_barcodes wrote no {prefix}.h5ad")
    return f
