"""Tests for flyte_bio.modules.gffread.

md5 anchors from the upstream gffread module snapshot (pinned version). Both
conversions are pure functions of input content, so both are md5-checked.
"""



from flyte_bio.modules.gffread import gffread_gff_to_gtf, gffread_transcripts_fasta
from tests.framework import assert_md5, env, fixture


@env.task
async def test_gff_to_gtf() -> None:
    # upstream case: gffread "sarscov2-gff3-gtf" (ext.args = -T)
    gff = await fixture("genomics/sarscov2/genome/genome.gff3")
    out = await gffread_gff_to_gtf(gff=gff)
    await assert_md5(out, "1ea0ae98d3388e0576407dc4a24ef428", label="gffread gff->gtf")


@env.task
async def test_transcripts_fasta() -> None:
    # upstream case: gffread "sarscov2-gff3-fasta" (ext.args = -w)
    gff = await fixture("genomics/sarscov2/genome/genome.gff3")
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    out = await gffread_transcripts_fasta(annotation=gff, fasta=fasta)
    await assert_md5(out, "5f8108fb51739a0588ccf0a251de919a", label="gffread transcripts fasta")


tests = [test_gff_to_gtf, test_transcripts_fasta]
