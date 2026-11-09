"""Tests for flyte_bio.modules.misopy.

The case and expected values are upstream's misopy/index test (modules
commit 9a1cc9f, misopy 0.5.4). As upstream, the index's files are sorted
into empty files, unstable ones (the shelves and pickles, checked by name)
and stable ones (``genes.gff`` and each chromosome directory, which upstream
hashes file by file: there, the gene's pickle).
"""

import asyncio
import hashlib
import os
import tarfile
import tempfile

from flyte_bio.modules.misopy import misopy_index
from tests.framework import env, fixture

EMPTY = [
    "compressed_ids_to_genes.shelve.bak",
    "compressed_ids_to_genes.shelve.dat",
    "compressed_ids_to_genes.shelve.dir",
]
PICKLE = "gene:ENSG00000198445.pickle"
UNSTABLE = [PICKLE, "genes_to_filenames.shelve.bak", "genes_to_filenames.shelve.dat", "genes_to_filenames.shelve.dir"]


def file_md5(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def classify_index(archive: str, prefix: str) -> tuple[list[str], list[str], dict[str, str], list[list[str]]]:
    """Unpack ``archive`` and sort ``<prefix>/``'s files as upstream's test: empty, unstable, stable files, dirs."""
    work = tempfile.mkdtemp()
    with tarfile.open(archive) as tar:
        tar.extractall(work, filter="data")
    index = os.path.join(work, prefix)
    assert os.path.isdir(index), f"the archive holds no {prefix}/"
    empty, unstable, stable_files, stable_dirs = [], [], {}, []
    for root, dirs, files in os.walk(index):
        for name in dirs:
            path = os.path.join(root, name)
            stable_dirs.append(sorted(f"{f}:md5,{file_md5(os.path.join(path, f))}" for f in os.listdir(path)))
        for name in files:
            path = os.path.join(root, name)
            if os.path.getsize(path) == 0:
                empty.append(name)
            elif "shelve" in name or name.endswith("pickle"):
                unstable.append(name)
            else:
                stable_files[name] = file_md5(path)
    return sorted(empty), sorted(unstable), stable_files, stable_dirs


@env.task
async def test_misopy_index() -> None:
    # upstream case: misopy/index "homo sapiens - gff3"
    archive = await misopy_index(await fixture("genomics/homo_sapiens/genome/genome.gff3"), prefix="homo_sapiens")
    assert archive.name == "homo_sapiens.tar.gz", f"unexpected name {archive.name}"
    empty, unstable, stable_files, stable_dirs = await asyncio.to_thread(
        classify_index, str(await archive.download()), "homo_sapiens"
    )
    assert empty == EMPTY, f"empty files {empty}"
    assert unstable == UNSTABLE, f"shelves and pickles {unstable}"
    assert stable_files == {"genes.gff": "6885b9a0b5bb2427a88805168b456d59"}, f"stable files {stable_files}"
    expected_dir = [f"{PICKLE}:md5,6f2b19f22d2265ab16a54d50916b5902"]
    assert stable_dirs == [expected_dir], f"chromosome directories {stable_dirs}"


tests = [test_misopy_index]
