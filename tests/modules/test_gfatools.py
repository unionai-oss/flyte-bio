"""Tests for flyte_bio.modules.gfatools.

The case and expected md5 are upstream's gfatools/gfa2fa test (modules commit
9a1cc9f, gfatools 0.5). nf-test computes the md5 over the decompressed FASTA.
"""

from flyte_bio.modules.gfatools import gfatools_gfa2fa
from tests.framework import assert_gunzipped_md5, env, fixture


@env.task
async def test_gfatools_gfa2fa() -> None:
    gfa = await fixture("genomics/sarscov2/illumina/gfa/assembly.gfa")
    fasta = await gfatools_gfa2fa(gfa, prefix="test")
    assert fasta.name == "test.fasta.gz", f"gfa2fa wrote {fasta.name}"
    await assert_gunzipped_md5(fasta, "ac1f38d05fdc1f0bba84f90daa854921", label="test.fasta.gz")


tests = [test_gfatools_gfa2fa]
