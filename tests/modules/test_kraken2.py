"""Tests for flyte_bio.modules.kraken2.

The cases and expected md5s are upstream's kraken2 tests (modules commit
9a1cc9f, kraken2 2.1.6). Upstream untars the test database in a setup step;
here :func:`fixture_dir` extracts it. Upstream snapshots only the report and
the per-read classifications; for the saved reads it checks their names, as
these tests do.
"""

import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.kraken2 import kraken2, kraken2_add, kraken2_build
from tests.framework import assert_md5, assert_nonempty, env, fixture, fixture_dir

FASTQ = "genomics/sarscov2/illumina/fastq/"
REPORT_MD5 = "4227755fe40478b8d7dc8634b489761e"


async def database() -> Dir:
    return await fixture_dir("genomics/sarscov2/genome/db/kraken2.tar.gz", "kraken2")


async def check_names(files: list[File], names: list[str], label: str) -> None:
    assert [f.name for f in files] == names, f"{label}: got {[f.name for f in files]}, expected {names}"
    for file in files:
        await assert_nonempty(file, label=f"{label} {file.name}")


@env.task
async def test_kraken2_single_end() -> None:
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await kraken2(reads, await database(), prefix="test", save_output_fastqs=True)
    await assert_md5(result.report, REPORT_MD5, label="report")
    await check_names(result.classified_reads_fastq, ["test.classified.fastq.gz"], "classified")
    await check_names(result.unclassified_reads_fastq, ["test.unclassified.fastq.gz"], "unclassified")
    assert result.classified_reads_assignment is None


@env.task
async def test_kraken2_paired_end() -> None:
    reads_1 = await fixture(FASTQ + "test_1.fastq.gz")
    reads_2 = await fixture(FASTQ + "test_2.fastq.gz")
    result = await kraken2(reads_1, await database(), reads_2=reads_2, prefix="test", save_output_fastqs=True)
    await assert_md5(result.report, REPORT_MD5, label="report")
    await check_names(
        result.classified_reads_fastq, ["test.classified_1.fastq.gz", "test.classified_2.fastq.gz"], "classified"
    )
    await check_names(
        result.unclassified_reads_fastq,
        ["test.unclassified_1.fastq.gz", "test.unclassified_2.fastq.gz"],
        "unclassified",
    )


@env.task
async def test_kraken2_single_end_save_reads_assignment() -> None:
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await kraken2(reads, await database(), prefix="test", save_reads_assignment=True)
    await assert_md5(result.report, REPORT_MD5, label="report")
    assert result.classified_reads_assignment is not None
    await assert_md5(result.classified_reads_assignment, "e7a90531f0d8d777316515c36fe4cae0", label="classifiedreads")
    assert result.classified_reads_fastq == [] and result.unclassified_reads_fastq == []
    assert await result.results.get_file("test.classified.fastq.gz") is None, "wrote reads without save_output_fastqs"


METAGENOME = "genomics/sarscov2/metagenome/"
TAXONOMY_MD5S = {
    "names.dmp": "130f9132095562e09c732679c562f5e9",
    "nodes.dmp": "c471c27a4ce85ae74d2c63633c9ce1e3",
    "prot.accession2taxid": "c0f96ba5dbb00150b4b805ba6dab7bea",
}


async def kraken2_test_library(custom_map: bool) -> File:
    """Upstream's kraken2/add test inputs, added to a library ``test.tar.gz``."""
    return await kraken2_add(
        [
            await fixture("genomics/sarscov2/genome/genome.fasta"),
            await fixture("genomics/sarscov2/genome/proteome.fasta"),
        ],
        taxonomy_names=await fixture(METAGENOME + "prot_names.dmp"),
        taxonomy_nodes=await fixture(METAGENOME + "prot_nodes.dmp"),
        accession2taxid=await gunzip(await fixture(METAGENOME + "prot.accession2taxid.gz")),
        seqid2taxid=await fixture(METAGENOME + "prot_seqid2taxid.map") if custom_map else None,
        prefix="test",
    )


async def kraken2_unpack(archive: File) -> Path:
    """Extract ``archive`` locally; returns the directory it was extracted into."""
    root = Path(tempfile.mkdtemp())
    with tarfile.open(await archive.download()) as tar:
        tar.extractall(root, filter="data")
    return root


def kraken2_md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


async def check_kraken2_add(custom_map: bool, map_md5: str | None) -> None:
    archive = await kraken2_test_library(custom_map)
    assert archive.name == "test.tar.gz", archive.name
    db = await kraken2_unpack(archive) / "test"
    added = sorted(p.name for p in (db / "library" / "added").iterdir())
    assert len(added) == 6, f"library/added: expected 6 files, got {added}"
    taxonomy = sorted(p.name for p in (db / "taxonomy").iterdir())
    assert taxonomy == sorted(TAXONOMY_MD5S), f"taxonomy files: {taxonomy}"
    for name, md5 in TAXONOMY_MD5S.items():
        actual = kraken2_md5(db / "taxonomy" / name)
        assert actual == md5, f"taxonomy/{name}: md5 {actual}, expected {md5}"
    if map_md5 is None:
        assert not (db / "seqid2taxid.map").exists(), "wrote a seqid2taxid.map without a custom one"
    else:
        actual = kraken2_md5(db / "seqid2taxid.map")
        assert actual == map_md5, f"seqid2taxid.map: md5 {actual}, expected {map_md5}"


@env.task
async def test_kraken2_add_protein_db() -> None:
    # upstream case: kraken2/add "sarscov2 protein_db"
    await check_kraken2_add(False, None)


@env.task
async def test_kraken2_add_protein_db_custom_seqid2taxid() -> None:
    # upstream case: kraken2/add "sarscov2 protein_db custom-seqid2map"
    await check_kraken2_add(True, "054c8ba74ccfca4412802066eac923f6")


async def check_kraken2_build(custom_map: bool, hash_md5: str, taxo_md5: str) -> None:
    result = await kraken2_build(await kraken2_test_library(custom_map), clean=True, args="--protein")
    names = sorted(f.name for f in await result.db.list_files())
    # Upstream names unmapped.txt in both snapshots, but nf-test's file() doesn't check that a file
    # exists: with the custom map every sequence maps, and kraken2 writes no unmapped.txt.
    expected = ["opts.k2d"] if custom_map else ["opts.k2d", "unmapped.txt"]
    for name in expected:
        assert name in names, f"kraken2_build wrote no {name}: {names}"
    hash_k2d = await result.db.get_file("hash.k2d")
    taxo_k2d = await result.db.get_file("taxo.k2d")
    assert hash_k2d is not None and taxo_k2d is not None, f"kraken2_build: {names}"
    await assert_md5(hash_k2d, hash_md5, label="hash.k2d")
    await assert_md5(taxo_k2d, taxo_md5, label="taxo.k2d")
    assert result.archive is None, "cleaned database came with an archive"


@env.task
async def test_kraken2_build_protein_db() -> None:
    # upstream case: kraken2/build "sarscov2 protein_db"
    await check_kraken2_build(False, "e9984a5e98f87c0488cb5e7618d5bbe0", "29d65b1796e09191fd7bdcaa24130459")


@env.task
async def test_kraken2_build_protein_db_custom_seqid2taxid() -> None:
    # upstream case: kraken2/build "sarscov2 protein_db - customs seqid2taxid"
    await check_kraken2_build(True, "d88e1e529cfdd8be4226a1181e766d4a", "24338e2d78f803f48bcc5653c6e51816")


tests = [
    test_kraken2_single_end,
    test_kraken2_paired_end,
    test_kraken2_single_end_save_reads_assignment,
    test_kraken2_add_protein_db,
    test_kraken2_add_protein_db_custom_seqid2taxid,
    test_kraken2_build_protein_db,
    test_kraken2_build_protein_db_custom_seqid2taxid,
]
