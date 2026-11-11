"""Tests for flyte_bio.modules.abacas.

The case is upstream's abacas test (modules commit 9a1cc9f, ABACAS 1.3.1)
with its config's ``-m -p nucmer``. As upstream, the nucmer deltas (which
record full paths) and the unused-contigs list are checked by name, the
other outputs by md5.
"""

from flyte_bio.modules.abacas import abacas
from tests.framework import assert_md5, env, fixture

EXPECTED_MD5S = {
    "test.abacas.MULTIFASTA.fa": "46c899ad70dcef8d14b5829fd8fbab82",
    "test.abacas.bin": "ec7f67eac462821df9b6912a908ba5ca",
    "test.abacas.crunch": "9a95358a9bd8ee97d1f2253d95623a17",
    "test.abacas.fasta": "5e6c403d3850d52f6bde956fa2403b13",
    "test.abacas.gaps": "5361af445b8d18a85c3af6527a97c89a",
    "test.abacas.gaps.tab": "193024ec9e5a553573519b218eb06e29",
    "test.abacas.nucmer.tiling": "0adaa0ce800d92c149a523d447148d95",
    "test.abacas.tab": "a5b9b452516f519a4201ff809655ef69",
}
NAMED_ONLY = ["test.abacas.nucmer.delta", "test.abacas.nucmer.filtered.delta", "test.abacas.unused.contigs.out"]


@env.task
async def test_abacas_scaffolds() -> None:
    # upstream case: abacas "sarscov2 - scaffolds_fasta - genome_fasta"
    result = await abacas(
        await fixture("genomics/sarscov2/illumina/fasta/scaffolds.fasta"),
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        prefix="test.abacas",
        args="-m -p nucmer",
    )
    names = sorted(f.name for f in await result.results.list_files())
    assert names == sorted([*EXPECTED_MD5S, *NAMED_ONLY]), f"abacas wrote {names}"
    for name, md5 in EXPECTED_MD5S.items():
        found = await result.results.get_file(name)
        assert found is not None, name
        await assert_md5(found, md5, label=name)


tests = [test_abacas_scaffolds]
