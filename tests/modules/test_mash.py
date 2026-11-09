"""Tests for flyte_bio.modules.mash.

The cases and expected md5s are upstream's mash/sketch, mash/dist and
mash/screen tests (modules commit 9a1cc9f, Mash 2.3). The screen case
screens the test reads against a sketch of the SARS-CoV-2 genome, made with
mash_sketch as in upstream's setup.
"""

from flyte_bio.modules.mash import mash_dist, mash_screen, mash_sketch
from tests.framework import assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
GENOME = "genomics/sarscov2/genome/genome.fasta"


@env.task
async def test_mash_sketch_fastq() -> None:
    # upstream case: mash/sketch "sarscov2 - illumina - fastq_gz"
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    sketch, stats = await mash_sketch(reads, prefix="test")
    await assert_md5(sketch, "fbd39f82fb7473ed33bb02ee04979a81", label="test.msh")
    await assert_md5(stats, "82532bf3f085c6c84b878e58d3b745a3", label="test.mash_stats")


@env.task
async def test_mash_dist() -> None:
    # upstream case: mash/dist "test-mash-dist"
    query = await fixture("genomics/sarscov2/illumina/fasta/contigs.fasta")
    dist = await mash_dist(query, await fixture(GENOME), prefix="test")
    assert dist.name == "test.txt", dist.name
    await assert_md5(dist, "1dfd5b53e1996f9fdb9d260b3e81c4d6", label="test.txt")


@env.task
async def test_mash_screen_fastq() -> None:
    # upstream case: mash/screen "sarscov2 - illumina - fastq_gz"
    sketch, _ = await mash_sketch(await fixture(GENOME), prefix="sars_db")
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    screen = await mash_screen(reads, sketch, prefix="test")
    assert screen.name == "test.screen", screen.name
    await assert_md5(screen, "ac8701e1aab651b2f36c6380b1351b11", label="test.screen")


tests = [test_mash_sketch_fastq, test_mash_dist, test_mash_screen_fastq]
