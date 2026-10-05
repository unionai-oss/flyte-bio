"""Tests for flyte_bio.modules.tx2gene.

Expected md5s are the upstream custom/tx2gene snapshot values; the output
is a pure function of the GTF and the transcript IDs in the quant file.
"""

from flyte_bio.modules.tx2gene import tx2gene
from tests.framework import assert_md5, env, fixture, fixture_dir

YEAST = "genomics/eukaryotes/saccharomyces_cerevisiae/"


@env.task
async def test_tx2gene() -> None:
    # upstream case: custom/tx2gene "saccharomyces_cerevisiae - gtf"
    gtf = await fixture(YEAST + "genome_gfp.gtf")
    quants = await fixture_dir(YEAST + "kallisto_results.tar.gz", "kallisto_results")
    table = await tx2gene(gtf, quants, quant_type="kallisto", gene_id="gene_id", extra="gene_name")
    await assert_md5(table, "0e2418a69d2eba45097ebffc2f700bfe", label="tx2gene")


@env.task
async def test_tx2gene_multiple_extra_attributes() -> None:
    # upstream case: custom/tx2gene "saccharomyces_cerevisiae - gtf - multiple extra attributes"
    gtf = await fixture(YEAST + "genome_gfp.gtf")
    quants = await fixture_dir(YEAST + "kallisto_results.tar.gz", "kallisto_results")
    table = await tx2gene(gtf, quants, quant_type="kallisto", gene_id="gene_id", extra="gene_name,gene_biotype")
    await assert_md5(table, "97223927dc2e0dae6c38bad96aaa6f49", label="tx2gene extra attributes")


tests = [test_tx2gene, test_tx2gene_multiple_extra_attributes]
