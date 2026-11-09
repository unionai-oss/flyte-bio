"""Tests for flyte_bio.modules.leafcutter.

The case and expected md5s are upstream's leafcutter/clusterregtools test
(modules commit 9a1cc9f, leafcutter 2.0.3): as upstream's setup, the
junctions of the chr6 BAM are extracted twice (as samples ``test1`` and
``test2``) with regtools, then clustered. nf-test hashes the ``.gz`` counts
decompressed.
"""

from flyte_bio.modules.leafcutter import leafcutter_clusterregtools
from flyte_bio.modules.regtools import regtools_junctionsextract
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.chr6.bam"


@env.task
async def test_leafcutter_clusterregtools() -> None:
    # upstream case: leafcutter/clusterregtools "homo_sapiens - bam - junction files"
    bam, bai = await fixture(BAM), await fixture(BAM + ".bai")
    junctions = [await regtools_junctionsextract(bam, bai, prefix=sample) for sample in ("test1", "test2")]
    result = await leafcutter_clusterregtools(junctions, prefix="test_cluster")
    await assert_md5(result.pooled, "ab5a2601994d8b34371fa49c36f7cf68", label="pooled")
    await assert_md5(result.refined, "46ad4cf57b90102241c317aea8b8cfc1", label="refined")
    await assert_md5(result.sortedlibs, "4014491357f92f44f8963a9d52c1e867", label="sortedlibs")
    await assert_gunzipped_md5(result.counts, "46e320419c32ecf890bd4909a1446ea0", label="perind counts")
    await assert_gunzipped_md5(result.numers, "3f276586d9d7431245392cbcbf96f826", label="perind numers")


tests = [test_leafcutter_clusterregtools]
