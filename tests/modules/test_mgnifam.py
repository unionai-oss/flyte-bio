"""Tests for flyte_bio.modules.mgnifam.

The cases and expected md5s are upstream's mgnifam/generatefamilies tests
(modules commit 9a1cc9f, mgnifam 4.0.0) on the proteinfamilies test data:
the plain run, the run with an SSI index from
:func:`flyte_bio.modules.hmmer.esl_sfetch_index` (as upstream's setup builds
it with hmmer/eslsfetchindex), and ``--prefetch_targets``. All three yield the
same nine families. nf-test hashes ``.gz`` outputs decompressed; upstream
checks only the log's name, as it records wall-clock timestamps.
"""

from flyte.io import File

from flyte_bio.modules.hmmer import esl_sfetch_index
from flyte_bio.modules.mgnifam import MgnifamFamilies, mgnifam_generatefamilies
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

DATA = "proteomics/proteinfamilies/"
CLUSTERING = DATA + "mgnify_proteins_trimmed_clustering.tsv"
FASTA = DATA + "mgnify_proteins.faa"

SEED_MSA = [
    "ca2b7d83b573f9ecea81bbadcf20056b",
    "efdb1a88939e93ccb27316f72fbf9877",
    "d4af6fb4656b796000c2e30435130977",
    "95b4dfeb2574190e99c1b3fd39625156",
    "c39515900a7fbfd3da4c565fbc8b10e1",
    "8f1ff690d44c3b98e7e6773455749664",
    "213351a99b1d69e5ac798fda5bb20623",
    "aa8ff1a005eeb2244b0ca6c5ee7edd9c",
    "e58ceca6683130fda471b0e4595869f1",
]
FULL_MSA = [
    "c5cbbba2a38f76a5c1cb68b07904adfe",
    "7d4d8120cf70ebf483fbb6dfa574b399",
    "714b712b687a38198e21ec38be27ee53",
    "6a23576bb5d34fb371e67c8fc7ef451a",
    "211f672d5cb6c30c7930e74757b126f0",
    "a53d033b230e1614625e15605270443d",
    "04a864e7c15c0dc4f3a07ef1cecf5827",
    "86f9d0f860a230f0adc2c0e7d04d1fc1",
    "e288d0d93e033dbd9e6c8fc04c9c7536",
]
HMM = [
    "678567b321cf7eee615cd379973ad4ee",
    "f43f95f769f88f7ee38a3ba9527842a8",
    "2d5b5f6ae86838ddb533cfe7f671b6cf",
    "adb16ddc40811ccf335701b1292eabfd",
    "238306dacf3f4c2bdc4b4713f28c2c4c",
    "1e9959d503bc262870558202ea3366e2",
    "149ca280f557496b02b345f825055668",
    "4b77146390d7a973c5caaceddf8a93e7",
    "3f1f5766c7674eb188147fb4eac31d4e",
]
RF = [
    "5c50c2b316dd413b3c3921f1bab84b8c",
    "a8de4ad7599e74d4e3459a1ca405cb1f",
    "bb6d3d37d8ad3aa983f7e85990991e0b",
    "197656079e2fc679fd987632d79dd2b8",
    "e766b28d59bc873f8e54d246bd375ce9",
    "dff2fdf7a9f45c94cd7c110487322c23",
    "36e272978769d4ca26f9bc545625021f",
    "be4c95323e3da4130ff938ae8e81efc3",
    "acdfbe7d162f80a5133e9b10dac869b2",
]
CHUNK_FILES = {
    "tsv": ("test_families.tsv", "7b8d5036a4cd176a6c28df840889e7ec"),
    "csv": ("test_metadata.csv", "ddf7014574136bf765c2b4df69a8926e"),
    "successful": ("test_successful.txt", "2a987074e3666639c71f1855313698c8"),
    "discarded": ("test_discarded.csv", "e22d5364cee2e59aa8d8dc1788e871e9"),
    "converged": ("test_converged.txt", "9c23343a9b58782cc1aab5046737db94"),
    "stats": ("test_mgnifam_stats.json", "e3bed1d11c7715f5a2204111dcddd82b"),
}


async def check_families(files: list[File], folder: str, suffix: str, md5s: list[str]) -> None:
    names = [f.name for f in files]
    expected = [f"test_{n}{suffix}" for n in range(1, len(md5s) + 1)]
    assert names == expected, f"{folder}: {names}, expected {expected}"
    for file, md5 in zip(files, md5s):
        check = assert_gunzipped_md5 if suffix.endswith(".gz") else assert_md5
        await check(file, md5, label=f"{folder}/{file.name}")


async def check(result: MgnifamFamilies) -> None:
    await check_families(result.seed_msa, "seed_msa", ".sto.gz", SEED_MSA)
    await check_families(result.full_msa, "full_msa", ".sto.gz", FULL_MSA)
    await check_families(result.hmm, "hmm", ".hmm.gz", HMM)
    await check_families(result.rf, "rf", ".txt", RF)
    for field, (name, md5) in CHUNK_FILES.items():
        file = getattr(result, field)
        assert file is not None and file.name == name, f"{field}: expected {name}, got {file}"
        await assert_md5(file, md5, label=name)
    assert result.reps_fasta is not None and result.reps_fasta.name == "test_reps.fasta.gz", result.reps_fasta
    await assert_gunzipped_md5(result.reps_fasta, "f9c52e85ecd9b7ca3651801584cebeb5", label="test_reps.fasta.gz")
    assert result.log is not None and result.log.name == "test.log", result.log


@env.task
async def test_mgnifam_generatefamilies() -> None:
    await check(await mgnifam_generatefamilies(await fixture(CLUSTERING), await fixture(FASTA), prefix="test"))


@env.task
async def test_mgnifam_generatefamilies_with_index() -> None:
    fasta = await fixture(FASTA)
    index = await esl_sfetch_index(fasta)
    await check(await mgnifam_generatefamilies(await fixture(CLUSTERING), fasta, fasta_index=index, prefix="test"))


@env.task
async def test_mgnifam_generatefamilies_prefetch_targets() -> None:
    result = await mgnifam_generatefamilies(
        await fixture(CLUSTERING), await fixture(FASTA), prefix="test", args="--prefetch_targets"
    )
    await check(result)


tests = [
    test_mgnifam_generatefamilies,
    test_mgnifam_generatefamilies_with_index,
    test_mgnifam_generatefamilies_prefetch_targets,
]
