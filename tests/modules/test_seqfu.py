"""Tests for flyte_bio.modules.seqfu.

The cases and expected md5s are upstream's seqfu/stats tests (modules commit
9a1cc9f, SeqFu 1.22.3).
"""

from flyte_bio.modules.seqfu import seqfu_stats
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/"


@env.task
async def test_seqfu_stats_faa() -> None:
    result = await seqfu_stats([await fixture(GENOME + "proteome.fasta")], prefix="test_faa")
    await assert_md5(result.stats, "0d6bf2cc788f7828761440a1689cac04", label="test_faa.tsv")
    await assert_md5(result.multiqc, "8f3c2edaf1ea5be912c9f99b21b2856c", label="test_faa_mqc.txt")


@env.task
async def test_seqfu_stats_multiple_files() -> None:
    files = [
        await fixture(GENOME + "proteome.fasta"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.gz"),
        await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz"),
    ]
    result = await seqfu_stats(files, prefix="test_multiple")
    await assert_md5(result.stats, "d016de3d84187a06c8e19b8dabccb3ae", label="test_multiple.tsv")
    await assert_md5(result.multiqc, "dbc6e762eebbf756cd0687807de60445", label="test_multiple_mqc.txt")


tests = [
    test_seqfu_stats_faa,
    test_seqfu_stats_multiple_files,
]
