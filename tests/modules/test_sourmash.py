"""Tests for flyte_bio.modules.sourmash.

The cases are upstream's sourmash/sketch tests (modules commit 9a1cc9f,
sourmash 4.9.4). Upstream doesn't snapshot the signatures; it checks only
that they record the staged file names, as these tests do.
"""

from flyte.io import File

from flyte_bio.modules.sourmash import sourmash_sketch
from tests.framework import env, fixture

SARSCOV2_FASTA = "genomics/sarscov2/genome/genome.fasta"
HAEMOPHILUS_FNA = "genomics/prokaryotes/haemophilus_influenzae/genome/genome.fna.gz"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(signatures: File, name: str, recorded: list[str]) -> None:
    assert signatures.name == name, f"got {signatures.name}, expected {name}"
    content = await text(signatures)
    for path in recorded:
        assert path in content, f"{name} doesn't record {path}"


@env.task
async def test_sourmash_single_fasta() -> None:
    signatures = await sourmash_sketch([await fixture(SARSCOV2_FASTA)], prefix="test")
    await check(signatures, "test.sig", ["library/genome.fasta"])


@env.task
async def test_sourmash_multiple_fasta_separate() -> None:
    library = [await fixture(SARSCOV2_FASTA), await fixture(HAEMOPHILUS_FNA)]
    signatures = await sourmash_sketch(library, prefix="test")
    await check(signatures, "test.sig", ["library/genome.fasta", "library/genome.fna.gz"])


@env.task
async def test_sourmash_multiple_fasta_merged() -> None:
    library = [await fixture(SARSCOV2_FASTA), await fixture(HAEMOPHILUS_FNA)]
    signatures = await sourmash_sketch(library, prefix="test", merge_sigs=True)
    await check(signatures, "test.sig", ["library/genome.f"])


@env.task
async def test_sourmash_merged_prefix_with_space() -> None:
    signatures = await sourmash_sketch([await fixture(SARSCOV2_FASTA)], prefix="test sample", merge_sigs=True)
    await check(signatures, "test sample.sig", ["library/genome.fasta"])


tests = [
    test_sourmash_single_fasta,
    test_sourmash_multiple_fasta_separate,
    test_sourmash_multiple_fasta_merged,
    test_sourmash_merged_prefix_with_space,
]
