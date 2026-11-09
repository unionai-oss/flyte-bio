"""Tests for flyte_bio.modules.carveme.

The cases and expected md5s are upstream's carveme/carve tests (modules
commit 9a1cc9f, CarveMe 1.6.6 with SCIP). As upstream, the model and the
gene-score and MILP-solution tables are checked for presence only; the MILP
problem and the protein and reaction scores by md5.
"""

from flyte_bio.modules.carveme import carveme_carve
from tests.framework import assert_md5, assert_nonempty, env, fixture


@env.task
async def test_carveme_carve_proteome() -> None:
    # upstream case: "candidatus_portiera_aleyrodidarum - proteome fasta - xml"
    fasta = await fixture("genomics/prokaryotes/candidatus_portiera_aleyrodidarum/genome/proteome.fasta")
    result = await carveme_carve(fasta, prefix="test")
    await assert_nonempty(result.model, label="test.xml")


@env.task
async def test_carveme_carve_mediadb_debug() -> None:
    # upstream case: "cazymes - faa - mediadb - debug - xml"
    fasta = await fixture("proteomics/database/protein_mini_with_cazymes.faa")
    mediadb = await fixture("generic/tsv/media_db.tsv")
    result = await carveme_carve(fasta, mediadb=mediadb, prefix="test", args="--debug --gapfill LB,M9")
    await assert_nonempty(result.model, label="test.xml")
    assert result.gene_scores is not None, "no test_gene_scores.tsv"
    assert result.milp_solution is not None, "no test_milp_solution.tsv"
    assert result.milp_problem is not None and result.protein_scores is not None and result.reaction_scores is not None
    await assert_md5(result.milp_problem, "1b8e947f09a05f4f8211192f92336e79", label="test_milp_problem.lp")
    await assert_md5(result.protein_scores, "e53ed545a67a62d853daf1629bc965c2", label="test_protein_scores.tsv")
    await assert_md5(result.reaction_scores, "d75d814c676619b77c1fcffa292c3e29", label="test_reaction_scores.tsv")


tests = [test_carveme_carve_proteome, test_carveme_carve_mediadb_debug]
