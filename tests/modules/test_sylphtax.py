"""Tests for flyte_bio.modules.sylphtax.

The case is upstream's sylphtax/taxprof test (modules commit 9a1cc9f,
sylph-tax 1.9.0), with its snapshot md5: SARS-CoV-2 reads profiled by sylph
against the SARS-CoV-2 genome, then classified with upstream's test taxonomy.
"""

from flyte_bio.modules.sylph import sylph_profile
from flyte_bio.modules.sylphtax import sylphtax_taxprof
from tests.framework import assert_md5, env, fixture

TAXONOMY = "https://github.com/nf-core/test-datasets/raw/taxprofiler/data/database/sylph/test_taxonomy.tsv.gz"


@env.task
async def test_sylphtax_taxprof_single_end() -> None:
    # upstream case: sylphtax/taxprof "sarscov2 illumina single-end [fastq_gz]"
    profile = await sylph_profile(
        await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz"),
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        prefix="test",
    )
    taxprof = await sylphtax_taxprof(profile, await fixture(TAXONOMY), prefix="test")
    assert taxprof.name == "test.sylphmpa", taxprof.name
    await assert_md5(taxprof, "1fd81c2f226e4766853d324d7aca1fc9", label="sylphmpa")


tests = [test_sylphtax_taxprof_single_end]
