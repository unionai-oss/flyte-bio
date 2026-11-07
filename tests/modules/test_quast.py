"""Tests for flyte_bio.modules.quast.

The cases and expected md5s are upstream's quast tests (modules commit
9a1cc9f, QUAST 5.3.0), including their counts of entries in the report
directory.
"""

import os

from flyte_bio.modules.quast import QuastResult, quast
from tests.framework import assert_md5, env, fixture

CONTIGS = "genomics/sarscov2/illumina/fasta/contigs.fasta"


async def check_report_entries(result: QuastResult, expected: int) -> None:
    local = await result.report.download()
    entries = sorted(os.listdir(local))
    assert len(entries) == expected, f"report has {len(entries)} entries, expected {expected}: {entries}"


@env.task
async def test_quast_reference() -> None:
    result = await quast(
        [await fixture(CONTIGS)],
        fasta=await fixture("genomics/sarscov2/genome/genome.fasta.gz"),
        gff=await fixture("genomics/sarscov2/genome/genome.gff3.gz"),
        prefix="test",
    )
    await check_report_entries(result, 15)
    await assert_md5(result.tsv, "1fc09ef5115692ed558cb31482bf3db0", label="test.tsv")
    assert result.misassemblies is not None and result.unaligned is not None
    await assert_md5(result.misassemblies, "b2c977af2f6e07603e0a0f4d1e32d058", label="test_misassemblies.tsv")
    await assert_md5(result.unaligned, "a3a36c1d206f89d12af74cef0c412f94", label="test_unaligned.tsv")


@env.task
async def test_quast_no_reference() -> None:
    result = await quast([await fixture(CONTIGS)], prefix="test")
    await check_report_entries(result, 12)
    await assert_md5(result.tsv, "6c4891c111a46ed20d5e3062b71aa390", label="test.tsv")
    assert result.misassemblies is None and result.unaligned is None and result.transcriptome is None


tests = [test_quast_reference, test_quast_no_reference]
