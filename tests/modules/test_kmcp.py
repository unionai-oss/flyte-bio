"""Tests for flyte_bio.modules.kmcp.

The cases and expected md5s are upstream's kmcp/compute and kmcp/index tests
(modules commit 9a1cc9f, KMCP 0.9.4), with their configs' prefixes
(``test_kmcp``). The compute "directory" case extracts upstream's test
archive with :func:`fixture_dir` (upstream untars it in a setup step). The
archives are extracted in the test and their files checked by name and md5.
Upstream also snapshots whether the index log contains "time elapsed"
(``false``), which isn't a property worth pinning, so it isn't ported.
"""

import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.kmcp import kmcp_compute, kmcp_index
from tests.framework import assert_md5, env, fixture, fixture_dir

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


tests = [test_kmcp_compute_directory, test_kmcp_compute_fasta, test_kmcp_index, test_kmcp_index_with_taxonomy]
