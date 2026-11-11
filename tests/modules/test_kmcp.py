"""Tests for flyte_bio.modules.kmcp.

The cases and expected md5s are upstream's kmcp/compute and kmcp/index tests
(modules commit 9a1cc9f, KMCP 0.9.4), with their configs' prefixes
(``test_kmcp``). The compute "directory" case extracts upstream's test
archive with :func:`fixture_dir` (upstream untars it in a setup step). The
archives are extracted in the test and their files checked by name and md5.
Upstream also snapshots whether the index log contains "time elapsed"
(``false``), which isn't a property worth pinning, so it isn't ported.

The search and profile cases are upstream's kmcp/search and kmcp/profile
tests (same commit), with their configs' prefixes and compute options
(``--by-seq`` for profile). The search results are gzipped, so their md5 is
over the decompressed file, as nf-test's. Upstream's "missing taxonomy" case
expects the run to fail with its error message; here the task fails with the
same message on stderr, which the error carries. The "wrong level" check is
upstream's stub case; it costs nothing here, so it is kept.
"""

import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.kmcp import kmcp_compute, kmcp_index, kmcp_profile, kmcp_search
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture, fixture_dir

GENOME = "genomics/sarscov2/genome/genome.fasta"
METAGENOME = "genomics/sarscov2/metagenome/"
INFO_MD5 = "484f52b6265e64bff2434d14698159bb"
COMPUTE_MD5S = {"_info.txt": INFO_MD5, "genome.fasta.unik": "158f49b1bdd8dd1aff91de66772509ff"}
INDEX_MD5S = {
    "R001/__db.yml": "17def8f599b7e841fd9851799ff052d1",
    "R001/__name_mapping.tsv": "82e6b631a4b9c5741cd75af717318ef9",
    "R001/_block001.uniki": "39000a457019ad6319667712e203a787",
}
TAXDMP_MD5S = {
    "taxdmp/genome2taxid.map": "43bef5503c1eccff6a746dea709a99e4",
    "taxdmp/names.dmp": "cac3fc251ff837dc31ef4b1832553d48",
    "taxdmp/nodes.dmp": "cf7b63f1194c1eb7baa5d3e1b745a69b",
}


async def check_kmcp_archive(archive: File, top: str, expected: dict[str, str]) -> None:
    """``archive`` holds ``top/`` with exactly the ``expected`` files (paths relative to it) and md5s."""
    assert archive.name == f"{top}.tar.gz", archive.name
    root = Path(tempfile.mkdtemp())
    with tarfile.open(await archive.download()) as tar:
        tar.extractall(root, filter="data")
    base = root / top
    found = {str(p.relative_to(base)): p for p in base.rglob("*") if p.is_file()}
    assert sorted(found) == sorted(expected), f"{archive.name}: {sorted(found)}"
    for name, md5 in expected.items():
        actual = hashlib.md5(found[name].read_bytes()).hexdigest()
        assert actual == md5, f"{name}: md5 {actual}, expected {md5}"


@env.task
async def test_kmcp_compute_directory() -> None:
    # upstream case: kmcp/compute "sarscov2 - directory"
    genomes = await fixture_dir("delete_me/kmcp/kmcp_profile.tar.gz", "kmcp_profile")
    result = await kmcp_compute(genomes, prefix="test_kmcp")
    await assert_md5(result.info, INFO_MD5, label="_info.txt")
    await check_kmcp_archive(result.archive, "test_kmcp", COMPUTE_MD5S)


@env.task
async def test_kmcp_compute_fasta() -> None:
    # upstream case: kmcp/compute "sarscov2 - fasta"
    result = await kmcp_compute(await fixture(GENOME), prefix="test_kmcp")
    await assert_md5(result.info, INFO_MD5, label="_info.txt")
    await check_kmcp_archive(result.archive, "test_kmcp", COMPUTE_MD5S)


@env.task
async def test_kmcp_index() -> None:
    # upstream case: kmcp/index "sarscov2 - fasta"
    compute = await kmcp_compute(await fixture(GENOME), prefix="test")
    result = await kmcp_index(compute.archive, prefix="test_kmcp")
    assert result.log.name == "test_kmcp.log", result.log.name
    await check_kmcp_archive(result.archive, "test_kmcp", INDEX_MD5S)


@env.task
async def test_kmcp_index_with_taxonomy() -> None:
    # upstream case: kmcp/index "sarscov2 - fasta - with taxonomy"
    compute = await kmcp_compute(await fixture(GENOME), prefix="test")
    seq2taxid = Path(tempfile.mkdtemp()) / "genome2taxid.map"
    seq2taxid.write_text("genome\t2697049\n")
    result = await kmcp_index(
        compute.archive,
        taxdmp=[await fixture(METAGENOME + "nodes.dmp"), await fixture(METAGENOME + "names.dmp")],
        seq2taxidmap=await File.from_local(str(seq2taxid)),
        prefix="test_kmcp",
    )
    await check_kmcp_archive(result.archive, "test_kmcp", INDEX_MD5S | TAXDMP_MD5S)


READS = "genomics/sarscov2/illumina/fastq/test_1.fastq.gz"


@env.task
async def test_kmcp_search() -> None:
    # upstream case: kmcp/search "sarscov2 - fasta"
    compute = await kmcp_compute(await fixture(GENOME), prefix="test")
    index = await kmcp_index(compute.archive, prefix="test_kmcp")
    result = await kmcp_search(await fixture(READS), index.archive, prefix="test")
    assert result.name == "test.gz", result.name
    await assert_gunzipped_md5(result, "84953d3d517f65722d43a2c3fdd04935", label="test.gz")


async def kmcp_profile_search(with_taxonomy: bool) -> tuple[File, File]:
    """Upstream's profile setup: the search results and the database archive."""
    compute = await kmcp_compute(await fixture(GENOME), prefix="test_kmcp", args="--by-seq")
    taxdmp, seq2taxid = None, None
    if with_taxonomy:
        taxdmp = [await fixture(METAGENOME + "nodes.dmp"), await fixture(METAGENOME + "names.dmp")]
        local = Path(tempfile.mkdtemp()) / "genome2taxid.map"
        local.write_text("MT192765.1\t2697049\n")
        seq2taxid = await File.from_local(str(local))
    index = await kmcp_index(compute.archive, taxdmp=taxdmp, seq2taxidmap=seq2taxid, prefix="test")
    return await kmcp_search(await fixture(READS), index.archive, prefix="test"), index.archive


@env.task
async def test_kmcp_profile_strain() -> None:
    # upstream case: kmcp/profile "sarscov2 - fasta - strain"
    search, db = await kmcp_profile_search(with_taxonomy=False)
    profile = await kmcp_profile(search, db, "strain", prefix="test_kmcp")
    assert profile.name == "test_kmcp.profile", profile.name
    await assert_md5(profile, "f5f32685a4da189844ded70d4dc9ef0b", label="test_kmcp.profile")


@env.task
async def test_kmcp_profile_species_missing_taxonomy() -> None:
    # upstream case: kmcp/profile "sarscov2 - fasta - species - missing taxonomy"
    search, db = await kmcp_profile_search(with_taxonomy=False)
    try:
        await kmcp_profile(search, db, "species", prefix="test_kmcp")
    except Exception as err:
        expected = "ERROR: --level 'species' requires taxdump files to be provided in the database directory."
        assert expected in str(err), f"kmcp profile failed without the taxonomy error: {err}"
    else:
        raise AssertionError("kmcp profile at species level succeeded without a taxonomy")


@env.task
async def test_kmcp_profile_species_with_taxonomy() -> None:
    # upstream case: kmcp/profile "sarscov2 - fasta - species - with taxonomy"
    search, db = await kmcp_profile_search(with_taxonomy=True)
    profile = await kmcp_profile(search, db, "species", prefix="test_kmcp")
    await assert_md5(profile, "edbc2d73fe97fb1bcf32beb0c4cca628", label="test_kmcp.profile")


@env.task
async def test_kmcp_profile_wrong_level() -> None:
    # upstream case (a stub there): kmcp/profile "sarscov2 - fasta - wrong level - stub"
    try:
        await kmcp_profile(await fixture(READS), await fixture(READS), "broken")
    except ValueError as err:
        assert "--level must be one of 'species', 'strain' or 'assembly" in str(err), str(err)
    else:
        raise AssertionError("kmcp_profile accepted level 'broken'")


tests = [
    test_kmcp_compute_directory,
    test_kmcp_compute_fasta,
    test_kmcp_index,
    test_kmcp_index_with_taxonomy,
    test_kmcp_search,
    test_kmcp_profile_strain,
    test_kmcp_profile_species_missing_taxonomy,
    test_kmcp_profile_species_with_taxonomy,
    test_kmcp_profile_wrong_level,
]
