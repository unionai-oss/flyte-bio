"""Tests for flyte_bio.modules.fusioninspector.

The case is upstream's fusioninspector test (modules commit 9a1cc9f,
FusionInspector 2.11.5): the starfusion test reads against upstream's
minigenome library (STAR-Fusion's build, see test_starfusion), inspecting
FusionInspector's ``test/fusion_targets.A.txt`` with ``--annotate
--examine_coding_effect``. Upstream checks the THRA--AC090627.1 evidence line
in both tables and snapshots the tables' and log's names and the contigs
GTF's md5; so do these tests.
"""

from flyte_bio.modules.fusioninspector import fusioninspector
from tests.framework import assert_md5, env, fixture
from tests.modules.test_starfusion import minigenome_lib

HUMAN = "genomics/homo_sapiens/genome/"
FUSION_TARGETS = "https://github.com/FusionInspector/FusionInspector/raw/master/test/fusion_targets.A.txt"
THRA_EVIDENCE = "THRA--AC090627.1\t72\t83\t72.00\t80.15"


@env.task
async def test_fusioninspector() -> None:
    # upstream case: fusioninspector "FUSIONINSPECTOR - test"
    result = await fusioninspector(
        reads_1=await fixture(HUMAN + "test_starfusion_rnaseq_1.fastq.gz"),
        reads_2=await fixture(HUMAN + "test_starfusion_rnaseq_2.fastq.gz"),
        fusion_list=await fixture(FUSION_TARGETS),
        genome_lib=await minigenome_lib(),
        prefix="test",
        args="--annotate --examine_coding_effect",
    )
    assert result.tsv is not None and result.tsv.name == "test.FusionInspector.fusions.tsv", result.tsv
    assert result.abridged_tsv is not None
    assert result.abridged_tsv.name == "test.FusionInspector.fusions.abridged.tsv", result.abridged_tsv.name
    for table in (result.tsv, result.abridged_tsv):
        async with table.open("rb") as fh:
            assert THRA_EVIDENCE in bytes(await fh.read()).decode(), f"{table.name} lacks {THRA_EVIDENCE!r}"
    assert result.log is not None and result.log.name == "FusionInspector.log", result.log
    assert [f.name for f in result.gtfs] == ["test.gtf"], [f.name for f in result.gtfs]
    await assert_md5(result.gtfs[0], "d523061bd8b443a014d4cca9406ec772", label="test.gtf")


tests = [test_fusioninspector]
