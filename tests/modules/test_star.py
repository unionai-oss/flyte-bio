"""Tests for flyte_bio.modules.star.

A STAR genome index isn't content-reproducible (the upstream snapshot only
md5s the stubbed/empty build), so this is a run-to-green check: the index
directory is produced and non-empty.
"""



from flyte_bio.modules.star import star_genome_generate
from tests.framework import assert_dir_nonempty, env, fixture


@env.task
async def test_genome_generate() -> None:
    # upstream case: star/genomegenerate "fasta_gtf"
    fasta = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    gtf = await fixture("genomics/homo_sapiens/genome/genome.gtf")
    index = await star_genome_generate(fasta=fasta, gtf=gtf)
    await assert_dir_nonempty(index, label="star genomeGenerate")


tests = [test_genome_generate]
