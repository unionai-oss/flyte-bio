"""Tests for flyte_bio.modules.medaka.

The case is upstream's medaka test (modules commit 9a1cc9f, medaka 2.2.2),
which snapshots the first 100 characters of the joined consensus lines.
"""

import gzip

from flyte_bio.modules.medaka import medaka
from tests.framework import env, fixture

CONSENSUS_START = ">MT192765.1GTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTC"


@env.task
async def test_medaka() -> None:
    reads = await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz")
    draft = await fixture("genomics/sarscov2/genome/genome.fasta")
    consensus = await medaka(reads, draft, prefix="test")
    assert consensus.path.endswith("test.fa.gz"), consensus.path
    async with consensus.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    start = "".join(lines)[:100]
    assert start == CONSENSUS_START, f"consensus starts {start!r}"


tests = [test_medaka]
