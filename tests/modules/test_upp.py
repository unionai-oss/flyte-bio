"""Tests for flyte_bio.modules.upp.

The cases are upstream's upp/align tests (modules commit 9a1cc9f, SEPP/UPP
4.5.6), with ``-m amino`` as upstream's test config. The first case's md5 is
upstream's snapshot; it relies on the module's default 4 CPUs, as upstream's
tests run with. The second extends a FAMSA backbone alignment of the seatoxin
sequences, with FAMSA's guide tree, by the toxin sequences; upstream treats
that alignment as unstable and checks only its name.
"""

from flyte_bio.modules.famsa import famsa_align, famsa_guidetree
from flyte_bio.modules.upp import upp_align
from tests.framework import assert_md5, assert_nonempty, env, fixture

MSA_DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/"
ARGS = "-m amino"


@env.task
async def test_upp_seatoxin_uncompressed() -> None:
    out = await upp_align(await fixture(MSA_DATA + "setoxin-ref.fa"), prefix="test", args=ARGS)
    assert out.name == "test.aln", f"output named {out.name}"
    await assert_md5(out, "e6b5291e9cdb40e9b7c72688e4da533b", label="test.aln")


@env.task
async def test_upp_with_tree_compressed() -> None:
    seatoxin = await fixture(MSA_DATA + "setoxin-ref.fa")
    tree = await famsa_guidetree(seatoxin, prefix="tree")
    backbone = await famsa_align(seatoxin, prefix="align_tree")
    out = await upp_align(
        await fixture(MSA_DATA + "toxin-ref.fa"),
        fasta_aligned=backbone,
        tree=tree,
        prefix="test_tree",
        compress=True,
        args=ARGS,
    )
    assert out.name == "test_tree.aln.gz", f"output named {out.name}"
    await assert_nonempty(out, label="test_tree.aln.gz")


tests = [test_upp_seatoxin_uncompressed, test_upp_with_tree_compressed]
