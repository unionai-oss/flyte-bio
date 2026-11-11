"""Tests for flyte_bio.modules.affy.

The cases and expected md5s are upstream's affy/justrma tests (modules commit
9a1cc9f): GSE38751's CEL files (upstream's untarred archive) and sample
sheet. As upstream, the RDS is checked by name only.
"""

from flyte_bio.modules.affy import affy_justrma
from tests.framework import assert_md5, env, fixture, fixture_dir

ARRAY = "genomics/homo_sapiens/array_expression/"
VERSIONS_MD5 = "d4eccefc5230f488e9204b46da35eaa2"


@env.task
async def test_affy_justrma_celfiles() -> None:
    # upstream case: "homo_sapiens genome affy_array [celfiles]"
    sheet = await fixture(ARRAY + "GSE38751.csv")
    celfiles = await fixture_dir(ARRAY + "GSE38751_RAW.tar")
    result = await affy_justrma(sheet, celfiles, prefix="test")
    await assert_md5(result.expression, "878dfdd16d266ba9b382a8c6f088908b", label="matrix")
    assert result.rds.name == "test_eset.rds", f"unexpected RDS {result.rds.name}"
    assert result.annotation is None, "annotation written without --build_annotation"
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


@env.task
async def test_affy_justrma_annotation() -> None:
    # upstream case: "homo_sapiens genome affy_array [celfiles] - annotation"
    sheet = await fixture(ARRAY + "GSE38751.csv")
    celfiles = await fixture_dir(ARRAY + "GSE38751_RAW.tar")
    result = await affy_justrma(
        sheet, celfiles, prefix="normalised", args="--sample_name_col name --build_annotation TRUE"
    )
    await assert_md5(result.expression, "462665542e9531f6d8993b745ce78ae4", label="matrix")
    assert result.annotation is not None and result.annotation.name == "hgu133a2.annotation.tsv"
    await assert_md5(result.annotation, "e2ade4cff2efaed4425a468110bbcd35", label="annotation")
    assert result.rds.name == "normalised_eset.rds", f"unexpected RDS {result.rds.name}"
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


tests = [test_affy_justrma_celfiles, test_affy_justrma_annotation]
