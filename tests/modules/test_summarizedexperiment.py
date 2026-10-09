"""Tests for flyte_bio.modules.summarizedexperiment.

The RDS embeds R/package serialization details, and the upstream snapshot
only checks part of the R session log, so this is a run-to-green check of
the upstream "unify" case: five gene-level assays plus row and column
metadata.
"""

from flyte_bio.modules.summarizedexperiment import summarized_experiment
from flyte_bio.modules.tx2gene import tx2gene
from flyte_bio.modules.tximport import tximport
from tests.framework import assert_nonempty, env, fixture, fixture_dir

YEAST = "genomics/eukaryotes/saccharomyces_cerevisiae/"


@env.task
async def test_summarized_experiment_unify() -> None:
    # upstream case: summarizedexperiment "multi_matrix - rowdata - coldata - unify summarized experiment"
    gtf = await fixture(YEAST + "genome_gfp.gtf")
    quants = await fixture_dir(YEAST + "kallisto_results.tar.gz", "kallisto_results")
    tx2gene_tsv = await tx2gene(gtf, quants, quant_type="kallisto", gene_id="gene_id", extra="gene_name")
    txi = await tximport(quants, tx2gene_tsv, quant_type="kallisto", prefix="test")
    samplesheet = await fixture(YEAST + "samplesheet.csv")
    rds = await summarized_experiment(
        {
            "counts": txi.counts_gene,
            "length_scale": txi.counts_gene_length_scaled,
            "counts_scaled": txi.counts_gene_scaled,
            "lengths": txi.lengths_gene,
            "tpm": txi.tpm_gene,
        },
        rowdata=tx2gene_tsv,
        coldata=samplesheet,
        prefix="gene",
    )
    await assert_nonempty(rds, label="SummarizedExperiment rds")


tests = [test_summarized_experiment_unify]
