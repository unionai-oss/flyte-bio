"""Tests for flyte_bio.modules.salmon.

A salmon index isn't content-reproducible, so this is a run-to-green check:
the index directory is produced and non-empty. Exercises the decoy-aware
path (genome FASTA supplied).
"""



from flyte_bio.modules.salmon import salmon_index
from tests.framework import assert_dir_nonempty, env, fixture


@env.task
async def test_index() -> None:
    # upstream case: salmon/index "sarscov2" (decoy-aware: genome + transcriptome)
    genome = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    transcripts = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    index = await salmon_index(transcript_fasta=transcripts, genome_fasta=genome)
    await assert_dir_nonempty(index, label="salmon index")


tests = [test_index]
