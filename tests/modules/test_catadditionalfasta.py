"""Tests for flyte_bio.modules.catadditionalfasta.

Both outputs are deterministic concatenations, so both are md5 anchors.
md5 from the upstream custom/catadditionalfasta snapshot (pinned version);
locally md5-verified against the vendored logic before wiring.
"""



from flyte_bio.modules.catadditionalfasta import cat_additional_fasta
from tests.framework import assert_md5, env, fixture


@env.task
async def test_cat_additional_fasta() -> None:
    # upstream case: custom/catadditionalfasta "sarscov2 - fastq - gtf"
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    gtf = await fixture("genomics/sarscov2/genome/genome.gtf")
    add = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    out_fasta, out_gtf = await cat_additional_fasta(fasta=fasta, gtf=gtf, add_fasta=add, biotype="test_biotype")
    await assert_md5(out_fasta, "6a20c1a2e465519320a0d01f338f5cb5", label="catadditionalfasta fasta")
    await assert_md5(out_gtf, "bc88d95e7f27540e6b9906105d5be361", label="catadditionalfasta gtf")


tests = [test_cat_additional_fasta]
