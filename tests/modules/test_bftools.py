"""Tests for flyte_bio.modules.bftools.

The case and expected md5 are upstream's bftools/showinf test (modules commit
9a1cc9f, Bio-Formats 8.0.0); nf-test hashes the OME-XML decompressed.
"""

from flyte_bio.modules.bftools import bftools_showinf
from tests.framework import assert_gunzipped_md5, env, fixture


@env.task
async def test_bftools_showinf_ome_tif() -> None:
    # upstream case: "ome tif"
    xml = await bftools_showinf(await fixture("imaging/ome-tiff/cycif-tonsil-cycle1.ome.tif"), prefix="0")
    assert xml.name == "0.xml.gz", f"got {xml.name}"
    await assert_gunzipped_md5(xml, "b4e08d86595e773c42483343b06cd51d", label="0.xml.gz")


tests = [test_bftools_showinf_ome_tif]
