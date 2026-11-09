"""Tests for flyte_bio.modules.miniasm.

The case and expected md5s are upstream's miniasm test (modules commit
9a1cc9f, miniasm 0.3-r179): the Bacteroides fragilis Nanopore reads and
their overlaps. Gzipped outputs are compared decompressed, as nf-test hashes
them.
"""

from flyte_bio.modules.miniasm import miniasm
from tests.framework import assert_gunzipped_md5, env, fixture

NANOPORE = "genomics/prokaryotes/bacteroides_fragilis/nanopore/"


@env.task
async def test_miniasm_bacteroides_fragilis() -> None:
    # upstream case: "bacteroides_fragilis - [fastq.gz, paf]- fasta"
    reads = await fixture(NANOPORE + "fastq/test.fastq.gz")
    paf = await fixture(NANOPORE + "overlap.paf")
    result = await miniasm(reads, paf, prefix="test")
    await assert_gunzipped_md5(result.gfa, "60e8b0c635905d20bc5296710749abde", label="gfa")
    await assert_gunzipped_md5(result.assembly, "3c0c2f48e6e366b43f7a12c43f7cb9e8", label="fasta")


tests = [test_miniasm_bacteroides_fragilis]
