"""Tests for flyte_bio.modules.orthofinder.

The cases and md5s are upstream's orthofinder tests (modules commit 9a1cc9f,
OrthoFinder 3.1.3). As upstream, every file in the results tree with one of
the stable names is md5-checked (the rest record timings or vary run to run),
and the species are named after their files (``sarscov2.fasta``,
``candidatus_portiera_aleyrodidarum.fasta``). The module's 4 CPUs match
upstream's test runs.
"""

import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.orthofinder import orthofinder
from tests.framework import assert_md5, env, fixture


def copy_as(local: str, name: str) -> str:
    target = Path(tempfile.mkdtemp()) / name
    target.write_bytes(Path(local).read_bytes())
    return str(target)


async def renamed(path: str, name: str) -> File:
    """The fixture at ``path``, staged under ``name`` (OrthoFinder names species by file)."""
    local = await (await fixture(path)).download()
    return await File.from_local(copy_as(str(local), name))


def extract_files(archive: str) -> list[Path]:
    root = Path(tempfile.mkdtemp())
    with tarfile.open(archive) as tar:
        tar.extractall(root, filter="data")
    return [p for p in root.rglob("*") if p.is_file()]


async def tar_gz_files(archive: File) -> list[Path]:
    return extract_files(str(await archive.download()))


async def check_stable_files(archive: File, expected: dict[str, str]) -> None:
    """Upstream's check: the md5 of every file in the results tree whose name is in ``expected``."""
    files = await tar_gz_files(archive)
    found = sorted((p.name, hashlib.md5(p.read_bytes()).hexdigest()) for p in files if p.name in expected)
    assert found == sorted(expected.items()), f"stable files: {found}"


@env.task
async def test_orthofinder_two_proteomes() -> None:
    # upstream case: orthofinder "sarscov2 - candidatus_portiera_aleyrodidarum - proteome"
    fastas = [
        await renamed("genomics/sarscov2/genome/proteome.fasta", "sarscov2.fasta"),
        await renamed(
            "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/genome/proteome.fasta",
            "candidatus_portiera_aleyrodidarum.fasta",
        ),
    ]
    result = await orthofinder(fastas, prefix="test")
    await check_stable_files(
        result.orthofinder,
        {
            "Statistics_PerSpecies.tsv": "984b5011a34d54527fe17896bfa36a2d",
            "SpeciesTree_Gene_Duplications_0.5_Support.txt": "8b7a673e2e8b6d1aeb697f2bb88afa18",
            "Orthogroups.tsv": "b06cb0585d8349f6138f46c753fe6fb5",
            "SpeciesTree_rooted.txt": "4d5ea525feebe479fca0c0768271ba81",
            "SpeciesTree_rooted_node_labels.txt": "2425ef0fe07d301411932eea4268060e",
        },
    )
    await assert_md5(result.orthogroups, "b06cb0585d8349f6138f46c753fe6fb5", label="Orthogroups.tsv")
    await assert_md5(result.sptree, "2425ef0fe07d301411932eea4268060e", label="SpeciesTree_rooted_node_labels.txt")
    # As upstream's snapshot: the single-copy directory exists but is empty.
    assert result.single_copy_seqs is not None, "no Single_Copy_Orthologue_Sequences"
    assert await tar_gz_files(result.single_copy_seqs) == [], "single-copy sequences found"


@env.task
async def test_orthofinder_resume() -> None:
    # upstream case: orthofinder "sarscov2 - candidatus_portiera_aleyrodidarum - proteome - resume"
    prior = await fixture("delete_me/orthofinder/WorkingDirectory.tar.gz")
    h1065 = await fixture(
        "https://raw.githubusercontent.com/nf-core/test-datasets/proteinfold/testdata/sequences/H1065.fasta"
    )
    result = await orthofinder([h1065], prior_run=prior, prefix="test_2")
    await check_stable_files(result.orthofinder, {"Orthogroups.tsv": "48cd9adbc60df119249251a3046b0ff8"})
    assert result.single_copy_seqs is None, "resume wrote single-copy sequences"


tests = [test_orthofinder_two_proteomes, test_orthofinder_resume]
