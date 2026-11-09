"""Tests for flyte_bio.modules.macsyfinder.

The cases are upstream's macsyfinder/download and macsyfinder/search tests
(modules commit 9a1cc9f, MacSyFinder 2.1.6), with their chain: the TXSScan
models are downloaded, then the proteome is searched, plain and gzipped
(with the test config's ``--db-type unordered``). As upstream, outputs are
checked by name: the search writes ``all_systems.tsv``, ``macsyfinder.out``
and ``macsyfinder.err`` and, for unordered replicons, no best solutions.
"""

import os
import tarfile

from flyte.io import File

from flyte_bio.modules.macsyfinder import MacsyfinderResult, macsyfinder_download, macsyfinder_search
from flyte_bio.modules.pigz import pigz_compress
from tests.framework import assert_nonempty, env, fixture

PROTEOME = "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/genome/proteome.fasta"


async def check_search(result: MacsyfinderResult) -> None:
    assert result.best_solutions == [], [os.path.basename(f.path) for f in result.best_solutions]
    for name, file in (
        ("all_systems.tsv", result.summary),
        ("macsyfinder.out", result.stdout),
        ("macsyfinder.err", result.stderr),
    ):
        assert file is not None, f"no {name}"
    assert result.summary is not None
    await assert_nonempty(result.summary, label="all_systems.tsv")


@env.task
async def test_macsyfinder_download() -> None:
    # upstream case: macsyfinder/download "macsyfinder_download - TXSScan"
    models = await macsyfinder_download("TXSScan")
    local = await models.download()
    with tarfile.open(local) as tar:
        names = tar.getnames()
    assert "models" in names, f"archive lacks models/: {names[:5]}"
    assert any(n.startswith("models/TXSScan/") for n in names), f"archive lacks models/TXSScan: {names[:5]}"


@env.task
async def test_macsyfinder_search_fasta() -> None:
    # upstream case: macsyfinder/search "macsyfinder_search - proteins - fasta"
    models = await macsyfinder_download("TXSScan")
    proteins: File = await fixture(PROTEOME)
    result = await macsyfinder_search(proteins, models, "TXSScan", prefix="test", args="--db-type unordered")
    await check_search(result)


@env.task
async def test_macsyfinder_search_fasta_gz() -> None:
    # upstream case: macsyfinder/search "macsyfinder_search - proteins - fasta.gz"
    models = await macsyfinder_download("TXSScan")
    proteins = await pigz_compress(await fixture(PROTEOME))
    result = await macsyfinder_search(proteins, models, "TXSScan", prefix="test", args="--db-type unordered")
    await check_search(result)


tests = [test_macsyfinder_download, test_macsyfinder_search_fasta, test_macsyfinder_search_fasta_gz]
