"""Tests for flyte_bio.modules.salmon.

salmon outputs aren't content-reproducible (the index, and quant results
embed run metadata), so these are run-to-green checks: the directories are
produced and non-empty. Exercises the decoy-aware index and reads-mode
quant; alignment-mode quant is covered via the pipeline sub-driver.
"""



from flyte_bio.modules.salmon import salmon_index, salmon_quant_reads
from tests.framework import assert_dir_nonempty, env, fixture


@env.task
async def test_index() -> None:
    # upstream case: salmon/index "sarscov2" (decoy-aware: genome + transcriptome)
    genome = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    transcripts = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    index = await salmon_index(transcript_fasta=transcripts, genome_fasta=genome)
    await assert_dir_nonempty(index, label="salmon index")


@env.task
async def test_quant_reads() -> None:
    # upstream case: salmon/quant "sarscov2 - single_end"
    genome = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    transcripts = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    index = await salmon_index(transcript_fasta=transcripts, genome_fasta=genome)
    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    gtf = await fixture("genomics/sarscov2/genome/genome.gtf")
    results = await salmon_quant_reads(reads=[reads], index=index, gtf=gtf)
    await assert_dir_nonempty(results, label="salmon quant results")


tests = [test_index, test_quant_reads]
