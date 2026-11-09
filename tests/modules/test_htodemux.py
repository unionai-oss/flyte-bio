"""Tests for flyte_bio.modules.htodemux.

The case and expected values are upstream's htodemux test (modules commit
9a1cc9f): upstream's Seurat object with an ``HTO`` assay. Every output is
compared by md5, the RDS included.
"""

from flyte_bio.modules.htodemux import htodemux
from tests.framework import assert_md5, env, fixture

SEURAT_OBJECT = "genomics/homo_sapiens/10xgenomics/cellranger/hashing_demultiplexing/htodemux.rds"


@env.task
async def test_htodemux_seurat_object() -> None:
    # upstream case: htodemux "seuratObject - rds"
    result = await htodemux(await fixture(SEURAT_OBJECT), prefix="test", assay="HTO")
    await assert_md5(result.assignment, "33d72b27c3d0b0f7394af18a18490a5c", label="assignment")
    await assert_md5(result.classification, "886f595f148d2ccdd5a9134b8eb815cb", label="classification")
    await assert_md5(result.params, "d01208cd0d0b76548cc71de04e76f3d2", label="params")
    await assert_md5(result.rds, "ec43d2d0217c3a87529f1eb4d0e05ffb", label="rds")
    await assert_md5(result.versions, "c34406eb377dbda4e62a36a1fb78350b", label="versions.yml")


tests = [test_htodemux_seurat_object]
