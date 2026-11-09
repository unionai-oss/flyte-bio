"""Tests for flyte_bio.modules.krakenuniq.

The cases and expected md5s are upstream's krakenuniq/build tests (modules
commit 9a1cc9f, KrakenUniq 1.0.4), with the test config's options (its
``--jellyfish-bin`` is passed by the module). Upstream builds over the NCBI
taxonomy that ``krakenuniq-download`` fetches; these tests use the test
data's SARS-CoV-2 protein taxonomy (``prot_names.dmp``, ``prot_nodes.dmp``:
the genome's lineage) instead, renamed ``names.dmp`` and ``nodes.dmp``. As
upstream, the build log, report and ``taxDB`` are checked by name only (they
hold timestamps, paths and the taxonomy's labels), every other file by md5.
"""

import shutil
import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.krakenuniq import krakenuniq_build
from tests.framework import assert_md5, env, fixture

METAGENOME = "genomics/sarscov2/metagenome/"
ARGS = "--work-on-disk --max-db-size 12 --kmer-len 15 --minimizer-len 13"
UNSTABLE = ["database-build.log", "database.report.tsv", "taxDB"]
STABLE_MD5S = {
    "database.idx": "f5298aa4215c956d2f737b5988721a24",
    "database.kdb": "46542be19979e867075f9e455a314e00",
    "database.kdb.counts": "a82149480a435210ec1e870f06b6cdb3",
    "database0.kdb": "46542be19979e867075f9e455a314e00",
}
INTERMEDIATE_MD5S = {
    "database.jdb": "48db5f6bb9b9e9db7f5bbdd0d8bec4a8",
    "database.kraken.tsv": "be9bd0cb3a8bea30989e65aa42a73727",
    "library-files.txt": "8a9817d2de8465eb436aa3aa8696a717",
}


async def krakenuniq_renamed(path: str, name: str) -> File:
    """The fixture at ``path``, uploaded again under ``name``."""
    local = Path(tempfile.mkdtemp()) / name
    shutil.copy(await (await fixture(path)).download(), local)
    return await File.from_local(str(local))


async def krakenuniq_test_db(keep_intermediate: bool) -> Dir:
    return await krakenuniq_build(
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        taxonomy=[
            await krakenuniq_renamed(METAGENOME + "prot_names.dmp", "names.dmp"),
            await krakenuniq_renamed(METAGENOME + "prot_nodes.dmp", "nodes.dmp"),
        ],
        seqid2taxid=await fixture(METAGENOME + "seqid2taxid.map"),
        keep_intermediate=keep_intermediate,
        prefix="test",
        args=ARGS,
    )


async def check_krakenuniq_db(db: Dir, stable: dict[str, str], unstable: list[str]) -> None:
    files = {f.name: f async for f in db.walk(recursive=True)}
    names = sorted(files)
    assert names == sorted([*stable, *unstable]), f"krakenuniq db: {names}"
    for name, md5 in stable.items():
        await assert_md5(files[name], md5, label=name)


@env.task
async def test_krakenuniq_build() -> None:
    # upstream case: krakenuniq/build "sarscov2 - fasta"
    await check_krakenuniq_db(await krakenuniq_test_db(False), STABLE_MD5S, ["taxDB"])


@env.task
async def test_krakenuniq_build_nocleanup() -> None:
    # upstream case: krakenuniq/build "sarscov2 - fasta - nocleanup"
    await check_krakenuniq_db(await krakenuniq_test_db(True), STABLE_MD5S | INTERMEDIATE_MD5S, UNSTABLE)


tests = [test_krakenuniq_build, test_krakenuniq_build_nocleanup]
