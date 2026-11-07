"""Tests for flyte_bio.modules.wget.

The case and expected md5 are upstream's wget test (modules commit 9a1cc9f,
wget 1.21.4): download a GFF3 from the test datasets, with upstream's test
prefix ``test``.
"""

from flyte_bio.modules.wget import wget
from tests.framework import assert_md5, env

GFF3 = "https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/modules/data/genomics/sarscov2/genome/genome.gff3"


@env.task
async def test_wget_sarscov2_gff() -> None:
    # upstream case: wget "sarscov2 - gff"
    gff = await wget(GFF3, suffix="gff3", prefix="test")
    assert gff.name == "test.gff3", f"downloaded {gff.name}"
    await assert_md5(gff, "357435a81a9981a0128e840ebe11051e")


tests = [
    test_wget_sarscov2_gff,
]
