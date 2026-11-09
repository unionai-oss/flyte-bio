"""Tests for flyte_bio.modules.kalign.

The cases and expected md5s are upstream's kalign/align tests (modules commit
9a1cc9f, Kalign 3.4.0); nf-test hashes ``.gz`` outputs decompressed. They rely
on the module's default 4 CPUs, as upstream's tests run with.
"""

from flyte_bio.modules.kalign import kalign_align
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

CONTIGS = "genomics/sarscov2/illumina/fasta/contigs.fasta"
CONTIGS_MD5 = "c165ecf48fb89862cc2a991cc3cadb2d"


@env.task
async def test_kalign_contigs_uncompressed() -> None:
    out = await kalign_align(await fixture(CONTIGS), prefix="test")
    assert out.name == "test.aln", f"output named {out.name}"
    await assert_md5(out, CONTIGS_MD5, label="test.aln")


@env.task
async def test_kalign_contigs_compressed() -> None:
    out = await kalign_align(await fixture(CONTIGS), prefix="test", compress=True)
    assert out.name == "test.aln.gz", f"output named {out.name}"
    await assert_gunzipped_md5(out, CONTIGS_MD5, label="test.aln.gz")


tests = [test_kalign_contigs_uncompressed, test_kalign_contigs_compressed]
