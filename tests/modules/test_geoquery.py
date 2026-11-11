"""Tests for flyte_bio.modules.geoquery.

The cases and expected md5s are upstream's geoquery/getgeo tests (modules
commit 9a1cc9f): series GSE50790, fetched from NCBI GEO. As upstream, the
RDS and the matrix are checked by name only (upstream marks the matrix
unstable).
"""

from flyte.io import File

from flyte_bio.modules.geoquery import geoquery_getgeo
from tests.framework import assert_md5, env

VERSIONS_MD5 = "6364f5393cdbde3b4c7c92aac6f70dd0"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_geoquery_getgeo() -> None:
    # upstream case: "test-geoquery-getgeo"
    result = await geoquery_getgeo("GSE50790")
    assert (result.rds.name, result.expression.name) == ("eset.rds", "matrix.tsv")
    await assert_md5(result.annotation, "4dbb874b49a11dff13557f58da643ee0", label="annotation")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")
    lines = [line.strip() for line in (await text(result.versions)).splitlines()]
    assert lines == ['"GEOQUERY_GETGEO":', "r-base: 4.3.2", "bioconductor-geoquery: 2.70.0"], lines


@env.task
async def test_geoquery_getgeo_with_metacols() -> None:
    # upstream case: "test-geoquery-getgeo-with-metacols"
    result = await geoquery_getgeo(
        "GSE50790", prefix="test.", args='--metacols "ID,ENTREZ_GENE_ID,Gene Symbol,Sequence Type"'
    )
    assert (result.rds.name, result.expression.name) == ("test.eset.rds", "test.matrix.tsv")
    await assert_md5(result.annotation, "40300ae222ae35fa54daf10a1b86e830", label="annotation")
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


tests = [test_geoquery_getgeo, test_geoquery_getgeo_with_metacols]
