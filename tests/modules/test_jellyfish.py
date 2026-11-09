"""Tests for flyte_bio.modules.jellyfish.

The cases are upstream's jellyfish/count tests (modules commit 9a1cc9f,
jellyfish 2.3.1). As upstream, the hash is unstable, so it is checked by
name (and here, that it is nonempty).
"""

from flyte.io import File

from flyte_bio.modules.jellyfish import jellyfish_count
from tests.framework import assert_nonempty, env, fixture


async def check(jf: File) -> None:
    assert jf.name == "test.jf", f"jellyfish wrote {jf.name}"
    await assert_nonempty(jf, label="test.jf")


@env.task
async def test_jellyfish_count_fasta() -> None:
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    await check(await jellyfish_count(fasta, kmer_length=7, size="100M", prefix="test"))


@env.task
async def test_jellyfish_count_fasta_gzipped() -> None:
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta.gz")
    await check(await jellyfish_count(fasta, kmer_length=7, size="100M", prefix="test"))


tests = [test_jellyfish_count_fasta, test_jellyfish_count_fasta_gzipped]
