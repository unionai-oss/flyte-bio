"""Tests for flyte_bio.modules.unicycler.

The case and expected md5s are upstream's unicycler test (modules commit
9a1cc9f, Unicycler 0.5.1): a long-read-only assembly of the Bacteroides
fragilis Nanopore reads. Gzipped outputs are compared decompressed, as
nf-test hashes them; as upstream, the log (which has timings) is only checked
for presence. The module's 4 CPUs match upstream's test CPUs.
"""

from flyte_bio.modules.unicycler import unicycler
from tests.framework import assert_gunzipped_md5, assert_nonempty, env, fixture


@env.task
async def test_unicycler_nanopore() -> None:
    # upstream case: "bacteroides_fragilis - nanopore"
    reads = await fixture("genomics/prokaryotes/bacteroides_fragilis/nanopore/fastq/test.fastq.gz")
    result = await unicycler(longreads=reads, prefix="test")
    await assert_gunzipped_md5(result.scaffolds, "8f2b98020361b047eec9dc4c23792d37", label="scaffolds.fa")
    await assert_gunzipped_md5(result.gfa, "849cec81b0ebf8b34554cf56ce0939d2", label="assembly.gfa")
    await assert_nonempty(result.log, label="unicycler.log")


tests = [test_unicycler_nanopore]
