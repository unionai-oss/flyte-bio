"""Tests for flyte_bio.modules.canu.

The case and expected md5s are upstream's canu test (modules commit 9a1cc9f,
Canu 2.3): HiCanu on the human HiFi test reads with upstream's test options,
which keep every stage to 1 thread and 7 GB. Upstream's test config caps the
task at 1 CPU, so the test pins 1 CPU, salting the cache (which ignores
resources). FASTAs are compared decompressed, as nf-test hashes them.

As in upstream's snapshot, Canu 2.3 writes no ``<prefix>.contigs.layout``,
and HiFi mode skips correction and trimming, so those outputs are absent.
"""

import flyte

from flyte_bio.modules.canu import canu_cmd, read_canu_results
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

ARGS = (
    "maxInputCoverage=800 minOverlapLength=1 hapThreads=1 corThreads=1 minInputCoverage=0 stopOnLowCoverage=0"
    " maxMemory=7 merylMemory=7 merylThreads=1 batMemory=7 batThreads=1 redMemory=7 redThreads=1 oeaMemory=7"
    " oeaThreads=1 obtovlThreads=1 utgovlThreads=1"
)


@env.task
async def test_canu_hicanu() -> None:
    # upstream case: "test_hicanu"
    reads = await fixture("genomics/homo_sapiens/pacbio/fastq/test_hifi.fastq.gz")
    task = canu_cmd.as_task().override(
        resources=flyte.Resources(cpu=1, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=1")
    )
    results = await task(reads=[reads], mode="-pacbio-hifi", genome_size="10k", prefix="test", args=ARGS)
    result = await read_canu_results(results, "test")
    assert result.assembly is not None, "canu wrote no contigs"
    await assert_gunzipped_md5(result.assembly, "aeb2768e013506b25b6778ddadb77466", label="contigs.fasta")
    await assert_gunzipped_md5(result.contigs, "8e108c4bcef54fef9f1cff3c8b53182e", label="unassembled.fasta")
    await assert_md5(result.report, "5108f1eaf984efa72b592507d6a5b1cb", label="report")
    assert result.contig_info is not None and result.contig_position is not None, "canu wrote no layout tables"
    await assert_md5(result.contig_info, "092286fa522a008e628b533eb56f7ff6", label="tigInfo")
    await assert_md5(result.contig_position, "ed2086b99f8f86607d8b1c793bae5096", label="readToTig")
    assert result.metadata is None, "canu wrote an unexpected contigs.layout"
    assert result.corrected_reads is None and result.corrected_trimmed_reads is None


tests = [test_canu_hicanu]
