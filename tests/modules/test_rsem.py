"""Tests for flyte_bio.modules.rsem.

The cases are upstream's rsem/preparereference and rsem/calculateexpression
tests (modules commit 9a1cc9f, RSEM 1.3.1, STAR 2.7.10a), with their snapshot
md5s: every file of the reference, and the gene and isoform results and the
model files of each expression estimate. As upstream's test config, the FASTQ
case builds the reference with ``--star`` and estimates with
``--star --star-gzipped-read-file``.
"""

from flyte_bio.modules.rsem import RsemExpression, rsem_calculateexpression, rsem_preparereference
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
TRANSCRIPTS_MD5 = "050c521a2719c2ae48267c1e65218f29"
INDEX_MD5S = {
    "genome.chrlist": "b190587cae0531f3cf25552d8aa674db",
    "genome.fasta": "f315020d899597c1b57e5fe9f60f4c3e",
    "genome.grp": "c2848a8b6d495956c11ec53efc1de67e",
    "genome.idx.fa": TRANSCRIPTS_MD5,
    "genome.n2g.idx.fa": TRANSCRIPTS_MD5,
    "genome.seq": "94da0c6b88c33e63c9a052a11f4f57c1",
    "genome.ti": "c9e4ae8d4d13a504eec2acf1b8589a66",
    "genome.transcripts.fa": TRANSCRIPTS_MD5,
}
GENES_MD5 = "c7ec226f76736ea805771e73553ae359"
ISOFORMS_MD5 = "99f7f80aa505b44ca429fdebbd7dd5d8"
STAT_MD5S = {
    "test.cnt": "76249e6b2f3c104f414aae596ba2c2f4",
    "test.model": "a7a4bc1734918ef5848604e3362b83e2",
    "test.theta": "de2e4490c98cc5383a86ae8225fd0a28",
}


@env.task
async def test_rsem_preparereference_homo_sapiens() -> None:
    # upstream case: rsem/preparereference "homo_sapiens"
    reference = await rsem_preparereference(
        await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.gtf")
    )
    files = {f.name: f for f in await reference.index.list_files()}
    assert sorted(files) == sorted(INDEX_MD5S), f"rsem reference holds {sorted(files)}"
    for name, md5 in INDEX_MD5S.items():
        await assert_md5(files[name], md5, label=name)
    assert reference.transcript_fasta.name == "genome.transcripts.fa", reference.transcript_fasta.name
    await assert_md5(reference.transcript_fasta, TRANSCRIPTS_MD5, label="transcript_fasta")


async def check_rsem_expression(expression: RsemExpression) -> None:
    await assert_md5(expression.counts_gene, GENES_MD5, label="genes.results")
    await assert_md5(expression.counts_transcript, ISOFORMS_MD5, label="isoforms.results")
    stat = {f.name: f for f in await expression.stat.list_files()}
    assert sorted(stat) == sorted(STAT_MD5S), f"rsem stat holds {sorted(stat)}"
    for name, md5 in STAT_MD5S.items():
        await assert_md5(stat[name], md5, label=name)


@env.task
async def test_rsem_calculateexpression_bam() -> None:
    # upstream case: rsem/calculateexpression "homo_sapiens - bam"
    reference = await rsem_preparereference(
        await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.gtf")
    )
    expression = await rsem_calculateexpression(
        await fixture("genomics/homo_sapiens/illumina/bam/rsem.transcript.bam"),
        reference.index,
        strandedness="forward",
        prefix="test",
    )
    await check_rsem_expression(expression)


@env.task
async def test_rsem_calculateexpression_fastq() -> None:
    # upstream case: rsem/calculateexpression "homo_sapiens - fastq"
    reference = await rsem_preparereference(
        await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.gtf"), args="--star"
    )
    expression = await rsem_calculateexpression(
        await fixture("genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz"),
        reference.index,
        reads_2=await fixture("genomics/homo_sapiens/illumina/fastq/test_rnaseq_2.fastq.gz"),
        strandedness="forward",
        prefix="test",
        args="--star --star-gzipped-read-file",
    )
    await check_rsem_expression(expression)
    assert expression.log is not None, "rsem-calculate-expression wrote no test.log"


tests = [
    test_rsem_preparereference_homo_sapiens,
    test_rsem_calculateexpression_bam,
    test_rsem_calculateexpression_fastq,
]
