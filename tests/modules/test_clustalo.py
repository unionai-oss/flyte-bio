"""Tests for flyte_bio.modules.clustalo.

The cases and expected md5s are upstream's clustalo/align tests (modules
commit 9a1cc9f, Clustal Omega 1.2.4). nf-test hashes ``.gz`` outputs
decompressed. The guide-tree case builds its tree with
:func:`clustalo_guidetree`, as upstream's does with clustalo/guidetree.
"""

from flyte.io import File

from flyte_bio.modules.clustalo import clustalo_align, clustalo_guidetree
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

CONTIGS = "genomics/sarscov2/illumina/fasta/contigs.fasta"
CONTIGS_MD5 = "74bb9a2820a91cf68db94dbd46787722"
PHYLOPLACE = "https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/phyloplace/testdata/"


def check_name(file: File, name: str) -> None:
    assert file.name == name, f"output named {file.name}, expected {name}"


@env.task
async def test_clustalo_contigs_uncompressed() -> None:
    out = await clustalo_align(await fixture(CONTIGS), prefix="test")
    check_name(out, "test.aln")
    await assert_md5(out, CONTIGS_MD5, label="test.aln")


@env.task
async def test_clustalo_contigs_compressed() -> None:
    out = await clustalo_align(await fixture(CONTIGS), prefix="test", compress=True)
    check_name(out, "test.aln.gz")
    await assert_gunzipped_md5(out, CONTIGS_MD5, label="test.aln.gz")


@env.task
async def test_clustalo_profile_hmm_in() -> None:
    sequences = await fixture(PHYLOPLACE + "PF14720_3_sequences.faa")
    out = await clustalo_align(sequences, hmm_in=await fixture(PHYLOPLACE + "PF14720.hmm"), prefix="test")
    await assert_md5(out, "9d3b05b9871b55ae5aaacedbc3f0f691", label="--hmm-in")


@env.task
async def test_clustalo_profile_profile1() -> None:
    sequences = await fixture(PHYLOPLACE + "PF14720_3_sequences.faa")
    out = await clustalo_align(sequences, profile1=await fixture(PHYLOPLACE + "PF14720_seed.alnfaa"), prefix="test")
    await assert_md5(out, "bc500eae88f3858f0f27800f8cd69698", label="--profile1")


@env.task
async def test_clustalo_contigs_guide_tree() -> None:
    contigs = await fixture(CONTIGS)
    tree = await clustalo_guidetree(contigs, prefix="test")
    out = await clustalo_align(contigs, tree=tree, prefix="test", compress=True)
    await assert_gunzipped_md5(out, CONTIGS_MD5, label="--guidetree-in")


tests = [
    test_clustalo_contigs_uncompressed,
    test_clustalo_contigs_compressed,
    test_clustalo_profile_hmm_in,
    test_clustalo_profile_profile1,
    test_clustalo_contigs_guide_tree,
]
