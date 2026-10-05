"""Tests for flyte_bio.modules.star.

STAR outputs aren't content-reproducible enough to md5 across thread counts
/ index parameters (the upstream snapshot md5s only the splice-junction
table), so these are run-to-green checks: the index and the alignment output
directory are produced and non-empty.
"""



from flyte_bio.modules.star import star_align, star_genome_generate
from tests.framework import assert_dir_nonempty, env, fixture

# Args from the upstream star/align "paired_end" test config.
ALIGN_ARGS = (
    "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate "
    "--outWigType bedGraph --outWigStrand Unstranded"
)


@env.task
async def test_genome_generate() -> None:
    # upstream case: star/genomegenerate "fasta_gtf"
    fasta = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    gtf = await fixture("genomics/homo_sapiens/genome/genome.gtf")
    index = await star_genome_generate(fasta=fasta, gtf=gtf)
    await assert_dir_nonempty(index, label="star genomeGenerate")


@env.task
async def test_align() -> None:
    # upstream case: star/align "homo_sapiens - paired_end"
    fasta = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    gtf = await fixture("genomics/homo_sapiens/genome/genome.gtf")
    index = await star_genome_generate(fasta=fasta, gtf=gtf, args="--genomeSAindexNbases 9")
    r1 = await fixture("genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz")
    r2 = await fixture("genomics/homo_sapiens/illumina/fastq/test_rnaseq_2.fastq.gz")
    out = await star_align(reads_1=r1, reads_2=[r2], index=index, gtf=gtf, args=ALIGN_ARGS)
    await assert_dir_nonempty(out, label="star align")


tests = [test_genome_generate, test_align]
