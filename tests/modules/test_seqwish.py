"""Tests for flyte_bio.modules.seqwish.

The cases and expected md5s are upstream's seqwish/induce tests (modules
commit 9a1cc9f, seqwish 0.7.11).
"""

from flyte_bio.modules.seqwish import seqwish_induce
from tests.framework import assert_md5, env, fixture


@env.task
async def test_seqwish_induce_transcriptome() -> None:
    # upstream case: seqwish/induce "seqwish_induce_transcriptome"
    paf = await fixture("genomics/sarscov2/genome/transcriptome.paf")
    fasta = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    result = await seqwish_induce([paf], fasta, prefix="test")
    assert result.gfa.path.endswith("/test.gfa"), result.gfa.path
    await assert_md5(result.gfa, "216a02d3aca322a457c31a62c628548e", label="test.gfa")


@env.task
async def test_seqwish_induce_pangenome() -> None:
    # upstream case: seqwish/induce "seqwish_induce_pangenome"
    paf = await fixture("pangenomics/homo_sapiens/pangenome.paf")
    fasta = await fixture("pangenomics/homo_sapiens/pangenome.fa")
    result = await seqwish_induce([paf], fasta, prefix="test")
    await assert_md5(result.gfa, "a73e998dc13a75d86d4685b12625d784", label="test.gfa")


tests = [test_seqwish_induce_transcriptome, test_seqwish_induce_pangenome]
