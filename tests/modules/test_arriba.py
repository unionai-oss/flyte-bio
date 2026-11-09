"""Tests for flyte_bio.modules.arriba.

The cases are upstream's arriba/arriba and arriba/visualisation tests
(modules commit 9a1cc9f, Arriba 2.5.0).

- arriba: upstream builds a STAR index of the homo_sapiens test genome
  (``--genomeSAindexNbases 9``), aligns ``test_rnaseq_{1,2}.fastq.gz`` with
  Arriba's recommended STAR options (chimeric reads ``WithinBAM``), and
  snapshots the md5s of both fusion tables. STAR runs through
  flyte_bio.modules.star with the test config's args and the read group
  upstream's star/align adds.
- visualisation: upstream snapshots only the version, so these check that a
  PDF is drawn: from the rnafusion test data's Arriba table, GTF, protein
  domains and cytobands; and from the arriba case's fusions with the test
  genome's GTF. For the latter, upstream takes the GRCh38 protein domains and
  cytobands from Arriba's 2.5.0 release archive (arriba/download, about
  500 MB); here they are the copies of those two files in the rnafusion test
  data.
"""

from flyte.io import File

from flyte_bio.modules.arriba import arriba_arriba, arriba_visualisation
from flyte_bio.modules.star import star_align, star_genome_generate
from tests.framework import assert_md5, env, fixture

HUMAN = "genomics/homo_sapiens/"
RNAFUSION = "https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/rnafusion/testdata/reference/"
# Upstream's test config for STAR_ALIGN, plus the read group upstream's module adds.
STAR_ARGS = (
    "--readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 "
    "--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 "
    "--alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip "
    "--chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 "
    "--chimSegmentReadGapMax 3 --chimMultimapNmax 50 --outSAMattrRGline ID:test SM:test"
)


async def arriba_test_fusions() -> tuple[File, File]:
    """Upstream's arriba case: STAR, then Arriba, on the homo_sapiens test data."""
    fasta = await fixture(HUMAN + "genome/genome.fasta")
    gtf = await fixture(HUMAN + "genome/genome.gtf")
    index = await star_genome_generate(fasta=fasta, gtf=gtf, args="--genomeSAindexNbases 9")
    aligned = await star_align(
        reads_1=await fixture(HUMAN + "illumina/fastq/test_rnaseq_1.fastq.gz"),
        reads_2=await fixture(HUMAN + "illumina/fastq/test_rnaseq_2.fastq.gz"),
        index=index,
        gtf=gtf,
        args=STAR_ARGS,
    )
    bam = await aligned.get_file("Aligned.out.bam")
    assert bam is not None, "STAR wrote no Aligned.out.bam"
    return await arriba_arriba(bam, fasta, gtf, prefix="test")


@env.task
async def test_arriba_arriba() -> None:
    # upstream case: arriba/arriba "arriba - homo_sapiens - paired_end"
    fusions, discarded = await arriba_test_fusions()
    assert fusions.name == "test.fusions.tsv", fusions.name
    assert discarded.name == "test.fusions.discarded.tsv", discarded.name
    await assert_md5(fusions, "7c3383f7eb6d79b84b0bd30a7ef02d70", label="fusions")
    await assert_md5(discarded, "445cb87a27f063e751f93498cf8d10b5", label="discarded fusions")


async def assert_pdf(pdf: File) -> None:
    assert pdf.name == "test.pdf", pdf.name
    async with pdf.open("rb") as fh:
        assert bytes(await fh.read(5)) == b"%PDF-", "not a PDF"


@env.task
async def test_arriba_visualisation() -> None:
    # upstream case: arriba/visualisation "ARRIBA_VISUALISATION"
    pdf = await arriba_visualisation(
        fusions=await fixture(RNAFUSION + "arriba/test_fastqs.arriba.fusions.tsv"),
        gtf=await fixture(RNAFUSION + "ensembl/Homo_sapiens.GRCh38.102.chr4.gtf"),
        prefix="test",
        protein_domains=await fixture(RNAFUSION + "arriba/protein_domains_hg38_GRCh38_v2.5.0.gff3"),
        cytobands=await fixture(RNAFUSION + "arriba/cytobands_hg38_GRCh38_v2.5.0.tsv"),
    )
    await assert_pdf(pdf)


@env.task
async def test_arriba_visualisation_arriba_test() -> None:
    # upstream case: arriba/visualisation "ARRIBA_VISUALISATION - arriba test"
    fusions, _ = await arriba_test_fusions()
    pdf = await arriba_visualisation(
        fusions=fusions,
        gtf=await fixture(HUMAN + "genome/genome.gtf"),
        prefix="test",
        protein_domains=await fixture(RNAFUSION + "arriba/protein_domains_hg38_GRCh38_v2.5.0.gff3"),
        cytobands=await fixture(RNAFUSION + "arriba/cytobands_hg38_GRCh38_v2.5.0.tsv"),
    )
    await assert_pdf(pdf)


tests = [test_arriba_arriba, test_arriba_visualisation, test_arriba_visualisation_arriba_test]
