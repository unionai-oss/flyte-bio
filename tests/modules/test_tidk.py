"""Tests for flyte_bio.modules.tidk.

The cases and md5s are upstream's tidk/explore, tidk/search and tidk/plot
tests (modules commit 9a1cc9f, tidk 0.2.7), with their test configs' ``args``
(explore: ``--minimum 5 --maximum 30``; search's bedgraph case:
``--extension bedgraph``).
"""

from flyte_bio.modules.tidk import tidk_explore, tidk_plot, tidk_search
from tests.framework import assert_md5, env, fixture

CHR21 = "genomics/homo_sapiens/genome/chr21/sequence/genome.fasta"
GENOME = "genomics/homo_sapiens/genome/genome.fasta"
EXPLORE_ARGS = "--minimum 5 --maximum 30"


@env.task
async def test_tidk_explore_chr21() -> None:
    # upstream case: tidk/explore "homo_sapiens-genome_21_fasta-success"
    result = await tidk_explore(await fixture(CHR21), prefix="test", args=EXPLORE_ARGS)
    await assert_md5(result.explore_tsv, "f7450552e49b5b58d5b85d268097a73b", label="test.tidk.explore.tsv")
    assert result.top_sequence is not None, "tidk explore found no top sequence"
    await assert_md5(result.top_sequence, "583dad19558da6e0e5f0b7e43840ed7c", label="test.top.sequence.txt")


@env.task
async def test_tidk_explore_no_top_sequence() -> None:
    # upstream case: tidk/explore "homo_sapiens-genome_fasta-no_top_sequence"
    result = await tidk_explore(await fixture(GENOME), prefix="test", args=EXPLORE_ARGS)
    await assert_md5(result.explore_tsv, "944f49381b397f9292fba77ac2fde706", label="test.tidk.explore.tsv")
    assert result.top_sequence is None, "tidk explore wrote a top sequence"


@env.task
async def test_tidk_search_bedgraph() -> None:
    # upstream case: tidk/search "homo_sapiens-genome_fasta-bedgraph"
    result = await tidk_search(await fixture(GENOME), "TTAGGG", prefix="test", args="--extension bedgraph")
    assert result.tsv is None, "tidk search wrote a TSV"
    assert result.bedgraph is not None, "tidk search wrote no bedgraph"
    await assert_md5(result.bedgraph, "8277cc74bf083f51584ba6d2b8a8013f", label="test.bedgraph")


@env.task
async def test_tidk_search_tsv() -> None:
    # upstream case: tidk/search "homo_sapiens-genome_fasta-tsv"
    result = await tidk_search(await fixture(GENOME), "TTAGGG", prefix="test")
    assert result.bedgraph is None, "tidk search wrote a bedgraph"
    assert result.tsv is not None, "tidk search wrote no TSV"
    await assert_md5(result.tsv, "5b44a89396f412eb571ea240ef9deedd", label="test.tsv")


@env.task
async def test_tidk_plot_chr21() -> None:
    # upstream case: tidk/plot "homo_sapiens-genome_21_fasta-success", from a search of chr21
    search = await tidk_search(await fixture(CHR21), "TTAGGG", prefix="test")
    assert search.tsv is not None, "tidk search wrote no TSV"
    svg = await tidk_plot(search.tsv, prefix="test")
    await assert_md5(svg, "d8f44a32f0b0dd9df235df15c779392c", label="test.svg")


tests = [
    test_tidk_explore_chr21,
    test_tidk_explore_no_top_sequence,
    test_tidk_search_bedgraph,
    test_tidk_search_tsv,
    test_tidk_plot_chr21,
]
