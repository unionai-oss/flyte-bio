"""Tests for flyte_bio.modules.cmaple.

The cases and expected md5s are upstream's cmaple tests (modules commit
9a1cc9f, CMAPLE 1.1.0): a protein alignment, and a gzipped DNA alignment
with a starting tree. Upstream snapshots the tree by md5 and only the log's
name, since the log records run times.
"""

from flyte_bio.modules.cmaple import CmapleResult, cmaple
from tests.framework import assert_md5, env, fixture

HAEMOPHILUS = "genomics/prokaryotes/haemophilus_influenzae/genome/"


async def check(result: CmapleResult, prefix: str, md5: str) -> None:
    assert result.treefile.name == f"{prefix}.treefile", result.treefile.name
    assert result.log.name == f"{prefix}.log", result.log.name
    await assert_md5(result.treefile, md5, label=f"{prefix}.treefile")


@env.task
async def test_cmaple_seatoxin_aa_aln() -> None:
    aln = await fixture("genomics/eukaryotes/anemonia_sulcata/seatoxin-ref.aln")
    await check(await cmaple(aln, prefix="test_aa"), "test_aa", "e8d45477e64ac8aae9ebda4859acdaf1")


@env.task
async def test_cmaple_haemophilus_dna_newick_aln_gz() -> None:
    aln = await fixture(HAEMOPHILUS + "genome.aln.gz")
    tree = await fixture(HAEMOPHILUS + "genome.aln.nwk")
    await check(await cmaple(aln, newick=tree, prefix="test_dna"), "test_dna", "bbc5f987e523bda0335ae7a2d70b982c")


tests = [test_cmaple_seatoxin_aa_aln, test_cmaple_haemophilus_dna_newick_aln_gz]
