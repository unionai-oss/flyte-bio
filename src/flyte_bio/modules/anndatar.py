"""anndatar — convert an AnnData h5ad to R objects.

Exposes :func:`anndatar_convert`, which reads an h5ad with the anndataR
package and saves it as a Seurat object (``<prefix>.seurat.rds``) and a
SingleCellExperiment (``<prefix>.sce.rds``), running upstream's script
(vendored) on upstream's image.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

ANNDATAR_IMAGE = (
    "community.wave.seqera.io/library/"
    "bioconductor-anndatar_bioconductor-rhdf5_bioconductor-singlecellexperiment_r-seurat:a0f51df063bb9b2a"
)

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

anndatar_convert_cmd = shell.create(
    name="anndatar_convert",
    image=ANNDATAR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "h5ad": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        S=({inputs.script}); H=({inputs.h5ad})
        Rscript "${S[0]}" "${H[0]}" {outputs.results}/{inputs.prefix}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "anndatar",
    anndatar_convert_cmd.as_task(),
)


@dataclass
class RdsObjects:
    seurat: File  # <prefix>.seurat.rds
    sce: File  # <prefix>.sce.rds


async def anndatar_convert(h5ad: File, prefix: str) -> RdsObjects:
    """``h5ad`` as ``<prefix>.seurat.rds`` and ``<prefix>.sce.rds``."""
    script = await File.from_local(str(path("anndatar_convert.r")))
    results = await anndatar_convert_cmd(script=script, h5ad=h5ad, prefix=prefix)

    async def pick(suffix: str) -> File:
        f = await results.get_file(f"{prefix}{suffix}")
        if f is None:
            raise FileNotFoundError(f"anndatar_convert wrote no {prefix}{suffix}")
        return f

    return RdsObjects(seurat=await pick(".seurat.rds"), sce=await pick(".sce.rds"))
