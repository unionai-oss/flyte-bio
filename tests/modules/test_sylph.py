"""Tests for flyte_bio.modules.sylph.

The cases are upstream's sylph/profile tests (modules commit 9a1cc9f, sylph
0.9.0): SARS-CoV-2 reads profiled against the SARS-CoV-2 genome. Upstream
snapshots only the table's header line, so the tests check that.
"""

from flyte.io import File

from flyte_bio.modules.sylph import sylph_profile
from tests.framework import env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
GENOME = "genomics/sarscov2/genome/genome.fasta"
HEADER = (
    "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov\t"
    "ANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\t"
    "Containment_ind\tNaive_ANI\tkmers_reassigned\tContig_name"
)


async def sylph_header(profile: File) -> str:
    assert profile.name == "test.tsv", profile.name
    async with profile.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()[0]


@env.task
async def test_sylph_profile_single_end() -> None:
    # upstream case: sylph/profile "sarscov2 illumina single-end [fastq_gz]"
    profile = await sylph_profile(await fixture(FASTQ + "test_1.fastq.gz"), await fixture(GENOME), prefix="test")
    header = await sylph_header(profile)
    assert header == HEADER, f"sylph profile header: {header!r}"


@env.task
async def test_sylph_profile_paired_end() -> None:
    # upstream case: sylph/profile "sarscov2 illumina paired-end [fastq_gz]"
    profile = await sylph_profile(
        await fixture(FASTQ + "test_1.fastq.gz"),
        await fixture(GENOME),
        reads_2=await fixture(FASTQ + "test_2.fastq.gz"),
        prefix="test",
    )
    header = await sylph_header(profile)
    assert header == HEADER, f"sylph profile header: {header!r}"


tests = [test_sylph_profile_single_end, test_sylph_profile_paired_end]
