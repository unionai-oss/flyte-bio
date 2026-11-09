"""Tests for flyte_bio.modules.traitar.

The cases are upstream's traitar/pfamget and traitar/run tests (modules
commit 9a1cc9f, Traitar 0.1.1), with their chain: Pfam is downloaded, then
phenotypes are predicted from a gzipped genome. Pfam-A (about 1.5 GiB) is
hashed streaming from disk; nf-test hashes the ``.gz`` decompressed, so both
files have upstream's one md5. As upstream, the predictions are checked by
name (upstream's snapshot of the single-vote table matches no output, so it
is empty there; here the table's presence is checked).
"""

import asyncio
import gzip
import hashlib
import os

from flyte.io import Dir, File

from flyte_bio.modules.traitar import traitar_pfamget, traitar_run
from tests.framework import env, fixture

PFAM_MD5 = "4bbbd0105827633e4933e344bfe0e92f"


def md5_of(path: str, unzip: bool) -> str:
    digest = hashlib.md5()
    with gzip.open(path, "rb") if unzip else open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


async def pfam_file(pfam: Dir, name: str) -> File:
    found = await pfam.get_file(name)
    assert found is not None, f"pfam_data lacks {name}"
    return found


@env.task
async def test_traitar_pfamget() -> None:
    # upstream case: traitar/pfamget "Download PFAM database - default version"
    pfam = await traitar_pfamget()
    names = sorted(os.path.basename(f.path) for f in await pfam.list_files())
    assert names == ["Pfam-A.hmm", "Pfam-A.hmm.gz"], names
    for name in names:
        local = str(await (await pfam_file(pfam, name)).download())
        actual = await asyncio.to_thread(md5_of, local, name.endswith(".gz"))
        assert actual == PFAM_MD5, f"{name}: md5 {actual}, expected {PFAM_MD5}"


@env.task
async def test_traitar_run_genome() -> None:
    # upstream case: traitar/run "traitar - genome - fna.gz"
    genome = await fixture("genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz")
    result = await traitar_run(genome, "from_nucleotides", await traitar_pfamget(), prefix="test_sample")
    assert os.path.basename(result.predictions_combined.path) == "predictions_majority-vote_combined.txt"
    assert os.path.basename(result.predictions_single_votes.path) == "predictions_single-votes_combined.txt"
    flat = sorted(os.path.basename(f.path) for f in result.predictions_flat)
    expected = ["predictions_flat_majority-votes_combined.txt", "predictions_flat_single-votes_combined.txt"]
    assert flat == expected, flat


tests = [test_traitar_pfamget, test_traitar_run_genome]
