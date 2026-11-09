"""Tests for flyte_bio.modules.sourmash.

The cases are upstream's sourmash/sketch, sourmash/index and sourmash/gather
tests (modules commit 9a1cc9f, sourmash 4.9.4). Upstream doesn't snapshot the
signatures; it checks only that they record the staged file names, as these
tests do. For the index it checks only the file name, as here (plus that it's
a nonempty file). Gather's CSV is checked by its decompressed md5, with
upstream's test options ``--threshold-bp 10 --ksize 31``.

The sourmash/compare case is upstream's at the same commit, with its test
config's ``--ksize 31``: two signatures of the same genome, compared into a
NumPy matrix, its labels and a CSV, checked by md5 as upstream.
"""

from flyte.io import File

from flyte_bio.modules.sourmash import sourmash_compare, sourmash_gather, sourmash_index, sourmash_sketch
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, env, fixture

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


@env.task
async def test_sourmash_index_sarscov2_genome() -> None:
    # upstream case: sourmash/index "sarscov2 genome [fasta]"
    signatures = await sourmash_sketch([await fixture(SARSCOV2_FASTA)], prefix="test")
    index = await sourmash_index([signatures], ksize=31, prefix="test")
    assert index.name == "test.sbt.zip", f"wrote {index.name}"
    await assert_nonempty(index, label="test.sbt.zip")


@env.task
async def test_sourmash_gather_sarscov2_genome() -> None:
    # upstream case: sourmash/gather "sarscov2 genome - fasta"
    query = await sourmash_sketch([await fixture(SARSCOV2_FASTA)], prefix="test1")
    database = await sourmash_sketch([await fixture(SARSCOV2_FASTA)], prefix="test2")
    result = await sourmash_gather(query, [database], prefix="test1", args="--threshold-bp 10 --ksize 31")
    assert result.result is not None, "gather wrote no test1.csv.gz"
    assert result.result.name == "test1.csv.gz", f"wrote {result.result.name}"
    await assert_gunzipped_md5(result.result, "969c138d20b84610a4f8a4867215b2d1", label="test1.csv.gz")
    unsaved = [result.unassigned, result.matches, result.prefetch, result.prefetch_csv]
    assert unsaved == [None] * 4, f"gather saved outputs it wasn't asked for: {unsaved}"


@env.task
async def test_sourmash_compare_sarscov2_genome() -> None:
    # upstream case: sourmash/compare "sarscov2 genome [fasta]"
    signatures = [
        await sourmash_sketch([await fixture(SARSCOV2_FASTA)], prefix="test1"),
        await sourmash_sketch([await fixture(SARSCOV2_FASTA)], prefix="test2"),
    ]
    result = await sourmash_compare(
        signatures, save_numpy_matrix=True, save_csv=True, prefix="group1", args="--ksize 31"
    )
    for file, name, md5 in [
        (result.matrix, "group1_comp.npy", "8a8683e1201544f28e2c2e21d2b89e4f"),
        (result.labels, "group1_comp.npy.labels.txt", "cd4807ea8634f741113647e656b70090"),
        (result.csv, "group1_comp.csv", "48366e16838afecfde44bff57ae9dddd"),
    ]:
        assert file is not None, f"compare wrote no {name}"
        assert file.name == name, f"wrote {file.name}, expected {name}"
        await assert_md5(file, md5, label=name)


tests = [
    test_sourmash_single_fasta,
    test_sourmash_multiple_fasta_separate,
    test_sourmash_multiple_fasta_merged,
    test_sourmash_merged_prefix_with_space,
    test_sourmash_index_sarscov2_genome,
    test_sourmash_gather_sarscov2_genome,
    test_sourmash_compare_sarscov2_genome,
]
