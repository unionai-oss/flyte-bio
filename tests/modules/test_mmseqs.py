"""Tests for flyte_bio.modules.mmseqs.

The cases and expected md5s are upstream's mmseqs/createdb, mmseqs/cluster,
mmseqs/linclust, mmseqs/createtsv, mmseqs/createindex and
mmseqs/tsv2exprofiledb tests (modules commit 9a1cc9f, MMseqs2 18.8cc5c). The
clustering tests run on 2 CPUs, as upstream's did: a clustering result has one
data file per thread, and upstream lists ``.0`` and ``.1``.
"""

import asyncio
import hashlib
import shutil
import tarfile
import tempfile
from pathlib import Path

import flyte
from flyte.io import Dir

from flyte_bio.modules.mmseqs import (
    mmseqs_cluster_cmd,
    mmseqs_createdb,
    mmseqs_createindex,
    mmseqs_createtsv,
    mmseqs_databases,
    mmseqs_linclust_cmd,
    mmseqs_taxonomy,
    mmseqs_tsv2exprofiledb,
)
from tests.framework import assert_md5, env, fixture

SARSCOV2_DB = "genomics/sarscov2/genome/db/mmseqs.tar.gz"
BACTEROIDES_GENOME = "genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz"
CLUSTER_PREFIX = "test_output_cluster"
CLUSTER_FILES = [f"{CLUSTER_PREFIX}.{suffix}" for suffix in ("0", "1", "dbtype", "index")]
TWO_CPUS = flyte.Resources(cpu=2, memory="4Gi")


async def dir_md5s(directory: Dir) -> dict[str, str]:
    """Each file's md5, by name."""
    md5s = {}
    async for file in directory.walk():
        async with file.open("rb") as fh:
            md5s[file.name] = hashlib.md5(bytes(await fh.read())).hexdigest()
    return md5s


async def check_dir(directory: Dir, names: list[str], md5s: dict[str, str], label: str) -> None:
    actual = await dir_md5s(directory)
    assert sorted(actual) == sorted(names), f"{label}: files {sorted(actual)}, expected {sorted(names)}"
    for name, md5 in md5s.items():
        assert actual[name] == md5, f"{label}: {name} md5 {actual[name]}, expected {md5}"


def extract_database(archive: str) -> Path:
    """The databases in ``archive``'s ``mmseqs/``, without its ``tmp/`` (which MMseqs2 doesn't read)."""
    root = Path(tempfile.mkdtemp())
    with tarfile.open(archive) as tar:
        tar.extractall(root, filter="data")
    shutil.rmtree(root / "mmseqs" / "tmp")
    return root / "mmseqs"


async def sarscov2_database() -> Dir:
    """Upstream's untarred ``mmseqs.tar.gz``: a sarscov2 database with its clustering and alignment results."""
    archive = await (await fixture(SARSCOV2_DB)).download()
    return await Dir.from_local(str(await asyncio.to_thread(extract_database, str(archive))))


CONTIGS_DB = {
    "test": "a2cda8768736a7a317a09d61556194bd",
    "test.dbtype": "4352d88a78aa39750bf70cd6f27bcaa5",
    "test.index": "4ba298b011e2472ce9f6b99fe6b6e3d5",
    "test.lookup": "32f88756dbcb6aaf7b239b0d61730f1b",
    "test.source": "9ada5b3ea6e1a7e16c4418eb98ae8d9d",
    "test_h": "21c399702a071bdeecce09f9d1df4531",
    "test_h.dbtype": "740bab4f9ec8808aedb68d6b1281aeb2",
    "test_h.index": "d767fb43b37c0a644c676b00f9f93477",
}
PROTEOME_DB = {
    "test": "1162504bc65aacf734abdcb0cdbe87de",
    "test.dbtype": "f1d3ff8443297732862df21dc4e57262",
    "test.index": "8cdcbc06c2b99fdb09f3d1735a76def9",
    "test.lookup": "3e27cb93d9ee875ad42a6f32f5651bdc",
    "test.source": "eaa64fc8a5f7ec1ee49b0dcbd1a72e9d",
    "test_h": "f258f8cc04f83c270a75e8b00a6d2d89",
    "test_h.dbtype": "740bab4f9ec8808aedb68d6b1281aeb2",
    "test_h.index": "844bf1950bcd37284fdc5d7117ee4241",
}


@env.task
async def test_mmseqs_createdb_contigs() -> None:
    db = await mmseqs_createdb(await fixture("genomics/sarscov2/illumina/fasta/contigs.fasta"), prefix="test")
    await check_dir(db, list(CONTIGS_DB), CONTIGS_DB, "createdb contigs")


@env.task
async def test_mmseqs_createdb_gzipped_proteome() -> None:
    db = await mmseqs_createdb(await fixture("genomics/sarscov2/genome/proteome.fasta.gz"), prefix="test")
    await check_dir(db, list(PROTEOME_DB), PROTEOME_DB, "createdb proteome.fasta.gz")


def two_cpu_task(cmd):
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    return cmd.as_task().override(resources=TWO_CPUS, cache=flyte.Cache(behavior="auto", salt="cpu=2"))


@env.task
async def test_mmseqs_cluster() -> None:
    task = two_cpu_task(mmseqs_cluster_cmd)
    result = await task(
        db=await sarscov2_database(), db_pattern="*.dbtype", prefix=CLUSTER_PREFIX, args="--remove-tmp-files 1 -v 3"
    )
    stable = {
        f"{CLUSTER_PREFIX}.dbtype": "b9d9c6dbc098c97ae446f612efd8eafd",
        f"{CLUSTER_PREFIX}.index": "9848b52b6df827d80a04f7c71c50056b",
    }
    await check_dir(result, CLUSTER_FILES, stable, "cluster")


@env.task
async def test_mmseqs_linclust() -> None:
    task = two_cpu_task(mmseqs_linclust_cmd)
    result = await task(
        db=await sarscov2_database(), db_pattern="*.dbtype", prefix=CLUSTER_PREFIX, args="--remove-tmp-files 1 -v 0"
    )
    stable = {
        f"{CLUSTER_PREFIX}.dbtype": "b9d9c6dbc098c97ae446f612efd8eafd",
        f"{CLUSTER_PREFIX}.index": "a0b78f31aee2c327d72f32919814baf1",
    }
    await check_dir(result, CLUSTER_FILES, stable, "linclust")


@env.task
async def test_mmseqs_createtsv_cluster() -> None:
    # As upstream, one Dir is the query, the target and (picked by its pattern) the clustering result.
    db = await sarscov2_database()
    tsv = await mmseqs_createtsv(db, db, db, prefix="test_query", result_pattern="*_clu.dbtype")
    await assert_md5(tsv, "c81449fb936b76aad6f925b965e84bc5", label="test_query.tsv")


@env.task
async def test_mmseqs_createtsv_taxonomy() -> None:
    query = await mmseqs_createdb(await fixture(BACTEROIDES_GENOME), prefix="test_query")
    silva = await mmseqs_databases("SILVA")
    taxonomy = await mmseqs_taxonomy(query, silva, prefix="test_query", args="--search-type 2")
    tsv = await mmseqs_createtsv(taxonomy, None, query, prefix="test_query")
    await assert_md5(tsv, "7e42aeb7daeeb54739e84fc6ecb0b5fc", label="test_query.tsv")


# Upstream snapshots the indexed directory's files but for these name endings:
# the index data itself (mmseqs.idx) and the per-thread files of the
# clustering and alignment results (only checked to exist).
CREATEINDEX_UNCHECKED = ["mmseqs.idx", *(f"_aln.{i}" for i in range(1, 8)), *(f"_clu.{i}" for i in range(1, 8))]
CREATEINDEX_FILES = {
    "mmseqs": "a5ba5a88d2995da3522dddb4abca07ef",
    "mmseqs.dbtype": "f1d3ff8443297732862df21dc4e57262",
    "mmseqs.fasta": "b40600ad3be77f076df716e6cf99c64c",
    "mmseqs.idx.dbtype": "9a198d4f48144e20661df7fd2dc41bf7",
    "mmseqs.idx.index": "3451d87bd8a5a7182d36d53edd213c59",
    "mmseqs.index": "c012bdab1c61eeafcb99d1b26650f3d0",
    "mmseqs.lookup": "fa898551a6b303614ae6e29c237b7fc6",
    "mmseqs.source": "16bef02c30aadbfa8d035596502f0aa2",
    "mmseqs.tsv": "4e7ba50ce2879660dc6595286bf0d097",
    "mmseqs_aln.0": "043e7e17756de6492a673c4b8c164892",
    "mmseqs_aln.dbtype": "cbe23d9655fccd4311acdd6261c81bd8",
    "mmseqs_aln.index": "9e17db59a56719569d20f72f0f48924d",
    "mmseqs_aln.tsv": "e6fc24b66126dec818f89491edb79fec",
    "mmseqs_clu.0": "86256b3bd5cd46a4ce865c477a381e08",
    "mmseqs_clu.dbtype": "b9d9c6dbc098c97ae446f612efd8eafd",
    "mmseqs_clu.index": "a0b78f31aee2c327d72f32919814baf1",
    "mmseqs_h": "2b2cb13de706b3bb48a38942d851f452",
    "mmseqs_h.dbtype": "740bab4f9ec8808aedb68d6b1281aeb2",
    "mmseqs_h.index": "0040e6a02964914a87ef1efbe9011cbf",
    "mmseqs_h.tsv": "8c733b4ef72b8c470ac30537cf6981bd",
    "mmseqs_profile": "6c29ffea085c8c78e97533f9bcbedd98",
    "mmseqs_profile.dbtype": "f2dd0dedb2c260419ece4a9e03b2e828",
    "mmseqs_profile.index": "e09d927ca8d4779b7cfb937f435cf7d3",
    "mmseqs_profile.lookup": "fa898551a6b303614ae6e29c237b7fc6",
    "mmseqs_profile.source": "16bef02c30aadbfa8d035596502f0aa2",
    "mmseqs_profile_h": "2b2cb13de706b3bb48a38942d851f452",
    "mmseqs_profile_h.dbtype": "740bab4f9ec8808aedb68d6b1281aeb2",
    "mmseqs_profile_h.index": "0040e6a02964914a87ef1efbe9011cbf",
    "mmseqs_seq": "e6c2c906a7ced527a75c7b175776f0b1",
    "mmseqs_seq.dbtype": "4352d88a78aa39750bf70cd6f27bcaa5",
    "mmseqs_seq.index": "09f8cf12406b814198156e56f5acd8dc",
    "mmseqs_seq.lookup": "fa898551a6b303614ae6e29c237b7fc6",
    "mmseqs_seq.source": "16bef02c30aadbfa8d035596502f0aa2",
    "mmseqs_seq.tsv": "fcc20038fa951bc085463bb49ffa5e25",
    "mmseqs_seq_h": "2b2cb13de706b3bb48a38942d851f452",
    "mmseqs_seq_h.dbtype": "740bab4f9ec8808aedb68d6b1281aeb2",
    "mmseqs_seq_h.index": "0040e6a02964914a87ef1efbe9011cbf",
}


@env.task
async def test_mmseqs_createindex_sarscov2() -> None:
    indexed = await mmseqs_createindex(await sarscov2_database(), args="--remove-tmp-files 1 -v 0")
    # The unchecked files are only listed: the index data runs to gigabytes even for this small database.
    files = [file async for file in indexed.walk()]
    unchecked = [file.name for file in files if any(file.name.endswith(end) for end in CREATEINDEX_UNCHECKED)]
    assert "mmseqs.idx" in unchecked, f"createindex: no mmseqs.idx in {sorted(file.name for file in files)}"
    checked = [file for file in files if file.name not in unchecked]
    assert sorted(file.name for file in checked) == sorted(CREATEINDEX_FILES), (
        f"createindex: files {sorted(file.name for file in checked)}"
    )
    for file in checked:
        await assert_md5(file, CREATEINDEX_FILES[file.name], label=f"createindex {file.name}")


# Upstream's listing also has the archive's tmp/, which the fixture drops.
EXPROFILE_NAMES = [
    *(f"mmseqs{suffix}" for suffix in ("", ".dbtype", ".fasta", ".index", ".lookup", ".source", ".tsv")),
    *(f"mmseqs_aln.{i}" for i in range(8)),
    "mmseqs_aln.dbtype",
    "mmseqs_aln.index",
    "mmseqs_aln.tsv",
    *(f"mmseqs_clu.{i}" for i in range(8)),
    "mmseqs_clu.dbtype",
    "mmseqs_clu.index",
    *(f"mmseqs_db{db}{suffix}" for db in ("", "_aln", "_h", "_seq", "_seq_h") for suffix in ("", ".dbtype", ".index")),
    *(f"mmseqs_h{suffix}" for suffix in ("", ".dbtype", ".index", ".tsv")),
    *(f"mmseqs_profile{suffix}" for suffix in ("", ".dbtype", ".index", ".lookup", ".source")),
    *(f"mmseqs_profile_h{suffix}" for suffix in ("", ".dbtype", ".index")),
    *(f"mmseqs_seq{suffix}" for suffix in ("", ".dbtype", ".index", ".lookup", ".source", ".tsv")),
    *(f"mmseqs_seq_h{suffix}" for suffix in ("", ".dbtype", ".index")),
]


@env.task
async def test_mmseqs_tsv2exprofiledb_sarscov2() -> None:
    # As upstream: many of the new databases are validly empty on this small
    # input, so the files are listed and a few checked by md5.
    result = await mmseqs_tsv2exprofiledb(await sarscov2_database())
    md5s = {
        "mmseqs.tsv": "4e7ba50ce2879660dc6595286bf0d097",
        "mmseqs": "a5ba5a88d2995da3522dddb4abca07ef",
        "mmseqs_seq.index": "09f8cf12406b814198156e56f5acd8dc",
    }
    await check_dir(result, EXPROFILE_NAMES, md5s, "tsv2exprofiledb")


tests = [
    test_mmseqs_createdb_contigs,
    test_mmseqs_createdb_gzipped_proteome,
    test_mmseqs_cluster,
    test_mmseqs_linclust,
    test_mmseqs_createtsv_cluster,
    test_mmseqs_createtsv_taxonomy,
    test_mmseqs_createindex_sarscov2,
    test_mmseqs_tsv2exprofiledb_sarscov2,
]
