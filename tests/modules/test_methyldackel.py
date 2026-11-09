"""Tests for flyte_bio.modules.methyldackel.

The cases and expected md5s are upstream's methyldackel/extract and
methyldackel/mbias tests (modules commit 9a1cc9f, MethylDackel 0.6.1), with
the extract test config's ``ext.args`` (none, then ``--methylKit``).
"""

from flyte.io import File

from flyte_bio.modules.methyldackel import methyldackel_extract, methyldackel_mbias
from tests.framework import assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/test.paired_end.methylated.sorted.bam"
GENOME = "genomics/sarscov2/genome/genome.fasta"


async def methyldackel_inputs() -> dict[str, File]:
    return {
        "bam": await fixture(BAM),
        "bai": await fixture(BAM + ".bai"),
        "fasta": await fixture(GENOME),
        "fai": await fixture(GENOME + ".fai"),
    }


@env.task
async def test_methyldackel_extract() -> None:
    # upstream case: "test-methyldackel-extract"
    result = await methyldackel_extract(**await methyldackel_inputs())
    assert result.methylkit == [], f"unexpected methylKit files {[f.name for f in result.methylkit]}"
    assert [f.name for f in result.bedgraph] == ["test.paired_end.methylated.sorted_CpG.bedGraph"], result.bedgraph
    await assert_md5(result.bedgraph[0], "abcc9d1db6e48d5aced567c337b563d4", label="bedGraph")


@env.task
async def test_methyldackel_extract_methylkit() -> None:
    # upstream case: "test-methyldackel-extract-methyl-kit"
    result = await methyldackel_extract(**await methyldackel_inputs(), args="--methylKit")
    assert result.bedgraph == [], f"unexpected bedGraph files {[f.name for f in result.bedgraph]}"
    assert [f.name for f in result.methylkit] == ["test.paired_end.methylated.sorted_CpG.methylKit"], result.methylkit
    await assert_md5(result.methylkit[0], "b1dc217909cd5f8174058c7af785f51c", label="methylKit")


@env.task
async def test_methyldackel_mbias() -> None:
    # upstream case: "test-methyldackel-mbias"
    txt = await methyldackel_mbias(**await methyldackel_inputs(), prefix="test")
    await assert_md5(txt, "357bb944dc2cdffcc47fa0d34376e016", label=txt.name)


tests = [test_methyldackel_extract, test_methyldackel_extract_methylkit, test_methyldackel_mbias]
