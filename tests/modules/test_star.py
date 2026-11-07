"""Tests for flyte_bio.modules.star.

STAR outputs aren't content-reproducible enough to md5 across thread counts
/ index parameters (the upstream snapshot md5s only the splice-junction
table), so these are run-to-green checks: the index and the alignment output
directory are produced and non-empty. STAR's index version is pure data, so
it is md5-checked.
"""

from flyte_bio.modules.star import star_align, star_genome_generate, star_indexversion
from tests.framework import assert_dir_nonempty, assert_md5, env, fixture

# Args from the upstream star/align "paired_end" test config.
ALIGN_ARGS = (
    "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outWigType bedGraph --outWigStrand Unstranded"
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
    out = await star_align(reads_1=r1, reads_2=r2, index=index, gtf=gtf, args=ALIGN_ARGS)
    await assert_dir_nonempty(out, label="star align")


@env.task
async def test_star_indexversion_normal_run() -> None:
    # upstream case: star/indexversion "normal run"
    version = await star_indexversion()
    assert version.name == "index_version.txt", version.name
    await assert_md5(version, "184a6b325b4ac98bd08ffc0bd8f046cd", label="index_version.txt")


tests = [test_genome_generate, test_align, test_star_indexversion_normal_run]
