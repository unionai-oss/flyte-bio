"""Tests for flyte_bio.modules.famsa.

The cases and expected md5s are upstream's famsa/align tests (modules commit
9a1cc9f, FAMSA 2.4.1); nf-test hashes ``.gz`` outputs decompressed. The
alignment depends on FAMSA's thread count, so these rely on the module's
default 4 CPUs, as upstream's tests run with. The guide-tree case builds its
tree with :func:`famsa_guidetree`, as upstream's does with famsa/guidetree.
"""

from flyte.io import File

from flyte_bio.modules.famsa import famsa_align, famsa_guidetree
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

CONTIGS = "genomics/sarscov2/illumina/fasta/contigs.fasta"
CONTIGS_MD5 = "561e958b5073a0adc976cea90173fc13"


def check_name(file: File, name: str) -> None:
    assert file.name == name, f"output named {file.name}, expected {name}"


@env.task
async def test_famsa_contigs_uncompressed() -> None:
    out = await famsa_align(await fixture(CONTIGS), prefix="test_uncompressed")
    check_name(out, "test_uncompressed.aln")
    await assert_md5(out, CONTIGS_MD5, label="test_uncompressed.aln")


@env.task
async def test_famsa_contigs_compressed() -> None:
    out = await famsa_align(await fixture(CONTIGS), prefix="test_compressed", compress=True)
    check_name(out, "test_compressed.aln.gz")
    await assert_gunzipped_md5(out, CONTIGS_MD5, label="test_compressed.aln.gz")


@env.task
async def test_famsa_contigs_guide_tree() -> None:
    contigs = await fixture(CONTIGS)
    tree = await famsa_guidetree(contigs, prefix="test_guide_tree")
    check_name(tree, "test_guide_tree.dnd")
    out = await famsa_align(contigs, tree=tree, prefix="test_guide_tree", compress=True)
    await assert_gunzipped_md5(out, CONTIGS_MD5, label="test_guide_tree.aln.gz")


tests = [test_famsa_contigs_uncompressed, test_famsa_contigs_compressed, test_famsa_contigs_guide_tree]
