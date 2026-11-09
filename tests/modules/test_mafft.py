"""Tests for flyte_bio.modules.mafft.

The cases and expected md5s are upstream's mafft/align tests (modules commit
9a1cc9f, MAFFT 7.520) and its mafft/guidetree test (MAFFT 7.525). nf-test
hashes ``.gz`` outputs decompressed.
"""

import asyncio
import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.mafft import mafft_align, mafft_guidetree
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

SCAFFOLDS = "genomics/sarscov2/illumina/fasta/scaffolds.fasta"
GENOME = "genomics/sarscov2/genome/genome.fasta"
INFORMATIVE_SITES = "genomics/sarscov2/genome/alignment/informative_sites.fas"
ALL_SITES = "genomics/sarscov2/genome/alignment/all_sites.fas"
CONTIGS = "genomics/sarscov2/illumina/fasta/contigs.fasta"
SCAFFOLDS_MD5 = "23426611f4a0df532b6708f072bd445b"


def gzip_copy(source: Path, name: str) -> Path:
    gzipped = Path(tempfile.mkdtemp()) / name
    gzipped.write_bytes(gzip.compress(source.read_bytes()))
    return gzipped


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


@env.task
async def test_mafft_scaffolds_uncompressed() -> None:
    out = await mafft_align(await fixture(SCAFFOLDS), prefix="test")
    assert file_name(out) == "test.fas"
    await assert_md5(out, SCAFFOLDS_MD5, label="test.fas")


@env.task
async def test_mafft_scaffolds_compressed() -> None:
    out = await mafft_align(await fixture(SCAFFOLDS), prefix="test", compress=True)
    assert file_name(out) == "test.fas.gz"
    await assert_gunzipped_md5(out, SCAFFOLDS_MD5, label="test.fas.gz")


@env.task
async def test_mafft_scaffolds_gzipped_input() -> None:
    # As upstream, the input is gzipped by the test; the alignment must equal the uncompressed run's.
    local = await (await fixture(SCAFFOLDS)).download()
    gzipped = await asyncio.to_thread(gzip_copy, Path(local), "scaffolds.fasta.gz")
    out = await mafft_align(await File.from_local(str(gzipped)), prefix="test")
    await assert_md5(out, SCAFFOLDS_MD5, label="test.fas")


@env.task
async def test_mafft_add_informative_sites() -> None:
    out = await mafft_align(await fixture(GENOME), add=await fixture(INFORMATIVE_SITES), prefix="test", compress=True)
    await assert_gunzipped_md5(out, "a57a34f1c566dea114dc1b13416536d4", label="--add")


@env.task
async def test_mafft_add_fragments_informative_sites() -> None:
    sites = await fixture(INFORMATIVE_SITES)
    out = await mafft_align(await fixture(GENOME), add_fragments=sites, prefix="test", compress=True)
    await assert_gunzipped_md5(out, "aed7f866c3a20dc9d2f2b4ad73515961", label="--addfragments")


@env.task
async def test_mafft_add_full_informative_sites() -> None:
    sites = await fixture(INFORMATIVE_SITES)
    out = await mafft_align(await fixture(GENOME), add_full=sites, prefix="test", compress=True)
    await assert_gunzipped_md5(out, "611cb0a65195a282f110f7f56e310c66", label="--addfull")


@env.task
async def test_mafft_add_profile_informative_sites() -> None:
    sites = await fixture(INFORMATIVE_SITES)
    out = await mafft_align(await fixture(GENOME), add_profile=sites, prefix="test", compress=True)
    await assert_gunzipped_md5(out, "c2b5caf39beff4473878e6aa4036ad43", label="--addprofile")


@env.task
async def test_mafft_add_long_informative_sites() -> None:
    sites = await fixture(INFORMATIVE_SITES)
    out = await mafft_align(await fixture(GENOME), add_long=sites, prefix="test", compress=True)
    await assert_gunzipped_md5(out, "e8868da70d1f3050a8daaee0e53b2fd9", label="--addlong")


@env.task
async def test_mafft_add_all_sites_fragments_informative_sites() -> None:
    out = await mafft_align(
        await fixture(GENOME),
        add=await fixture(ALL_SITES),
        add_fragments=await fixture(INFORMATIVE_SITES),
        prefix="test",
        compress=True,
    )
    await assert_gunzipped_md5(out, "aed7f866c3a20dc9d2f2b4ad73515961", label="--add --addfragments")


@env.task
async def test_mafft_guidetree_contigs() -> None:
    tree = await mafft_guidetree(await fixture(CONTIGS), prefix="test")
    assert file_name(tree) == "test.dnd", f"output named {file_name(tree)}"
    await assert_md5(tree, "e376b2b428c095b1d41c26bc0e3789ba", label="test.dnd")


tests = [
    test_mafft_scaffolds_uncompressed,
    test_mafft_scaffolds_compressed,
    test_mafft_scaffolds_gzipped_input,
    test_mafft_add_informative_sites,
    test_mafft_add_fragments_informative_sites,
    test_mafft_add_full_informative_sites,
    test_mafft_add_profile_informative_sites,
    test_mafft_add_long_informative_sites,
    test_mafft_add_all_sites_fragments_informative_sites,
    test_mafft_guidetree_contigs,
]
