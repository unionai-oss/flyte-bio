"""Tests for flyte_bio.modules.multiseqdemux.

The case and expected values are upstream's multiseqdemux test (modules
commit 9a1cc9f): upstream's Seurat object with an ``HTO`` assay. Every output
is compared by md5, the RDS included.
"""

from flyte_bio.modules.multiseqdemux import multiseqdemux
from tests.framework import assert_md5, env, fixture

SEURAT_OBJECT = "genomics/homo_sapiens/10xgenomics/cellranger/hashing_demultiplexing/htodemux.rds"


@env.task
async def test_multiseqdemux_seurat_object() -> None:
    # upstream case: multiseqdemux "seuratObject - rds"
    result = await multiseqdemux(await fixture(SEURAT_OBJECT), prefix="test", assay="HTO")
    await assert_md5(result.params, "03b37983ec6c4e9eabb789a81e936a4e", label="params")
    await assert_md5(result.results, "9d4475b8c2778e7e361747dd521072e4", label="results")
    await assert_md5(result.rds, "f9b34865c0845ffaa1e0b0447adda202", label="rds")
    await assert_md5(result.versions, "010faafd9d09f3a732eca24e43b9abf7", label="versions.yml")


tests = [test_multiseqdemux_seurat_object]
