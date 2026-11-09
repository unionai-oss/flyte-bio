"""Tests for flyte_bio.modules.tcoffee.

The cases are upstream's tcoffee/align, tcoffee/alncompare,
tcoffee/consensus, tcoffee/extractfrompdb, tcoffee/irmsd, tcoffee/regressive,
tcoffee/seqreformat and tcoffee/tcs tests (modules commit 9a1cc9f, T-Coffee
13.46.0.919e8c6b) on the seatoxin data, with the ``args`` and prefixes of
upstream's test configs. (The seqreformat and alncompare tests don't load
their ``nextflow.config``, so they run with no ``args``; for alncompare that is
the same as its config's ``-compare_mode tc``.) Where upstream snapshots an
output's md5 (alncompare, extractfrompdb, seqreformat) these assert it; where
it checks only that the output holds a sequence (``1ahl``), so do these.

Fixture steps follow upstream's setups: guide trees from
:func:`famsa_guidetree` (align, regressive) and :func:`clustalo_guidetree`
(consensus), the consensus's input alignments from :func:`famsa_align` and
:func:`magus_align`, the TCS library from :func:`tcoffee_align`. Upstream's
pigz/compress setup steps are replaced by gzipping in the test.
"""

import asyncio
import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.clustalo import clustalo_guidetree
from flyte_bio.modules.famsa import famsa_align, famsa_guidetree
from flyte_bio.modules.magus import magus_align
from flyte_bio.modules.tcoffee import (
    tcoffee_align,
    tcoffee_alncompare,
    tcoffee_consensus,
    tcoffee_extractfrompdb,
    tcoffee_irmsd,
    tcoffee_regressive,
    tcoffee_seqreformat,
    tcoffee_tcs,
)
from tests.framework import assert_md5, assert_nonempty, env, fixture, fixture_dir

MSA_DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/"
SEQUENCES = MSA_DATA + "setoxin-ref.fa"
REFERENCE_MSA = MSA_DATA + "setoxin.ref"
STRUCTURES = MSA_DATA + "structures/seatoxin-ref.tar.gz"
TEMPLATE = MSA_DATA + "templates/seatoxin-ref_template.txt"


def gzip_copy(source: Path, name: str) -> Path:
    gzipped = Path(tempfile.mkdtemp()) / name
    gzipped.write_bytes(gzip.compress(source.read_bytes()))
    return gzipped


async def gzipped(file: File) -> File:
    """``file`` gzipped, as ``<name>.gz``: what upstream's pigz/compress setup step makes."""
    local = await file.download()
    out = await asyncio.to_thread(gzip_copy, Path(local), file.name + ".gz")
    return await File.from_local(str(out))


async def text_of(file: File) -> str:
    async with file.open("rb") as fh:
        raw = bytes(await fh.read())
    return (gzip.decompress(raw) if file.name.endswith(".gz") else raw).decode()


async def check_contains(file: File, name: str, wanted: str = "1ahl") -> None:
    assert file.name == name, f"output named {file.name}, expected {name}"
    assert wanted in await text_of(file), f"{name} lacks {wanted}"


async def seatoxin_structures() -> list[File]:
    structures = await fixture_dir(STRUCTURES, "seatoxin-ref")
    return sorted([f async for f in structures.walk()], key=lambda f: f.name)


# align


@env.task
async def test_tcoffee_align_sequence_compressed() -> None:
    result = await tcoffee_align(await fixture(SEQUENCES), prefix="test", compress=True)
    await check_contains(result.alignment, "test.aln.gz")


@env.task
async def test_tcoffee_align_sequence_uncompressed() -> None:
    result = await tcoffee_align(await fixture(SEQUENCES), prefix="test")
    await check_contains(result.alignment, "test.aln")


@env.task
async def test_tcoffee_align_with_guide_tree() -> None:
    sequences = await fixture(SEQUENCES)
    tree = await famsa_guidetree(sequences, prefix="test")
    result = await tcoffee_align(sequences, tree=tree, prefix="test", compress=True)
    await check_contains(result.alignment, "test.aln.gz")


@env.task
async def test_tcoffee_align_with_structure() -> None:
    result = await tcoffee_align(
        await fixture(SEQUENCES),
        accessory=await seatoxin_structures(),
        prefix="test",
        compress=True,
        args="-method TMalign_pair -output fasta_aln",
    )
    await check_contains(result.alignment, "test.aln.gz")


@env.task
async def test_tcoffee_align_with_lib() -> None:
    result = await tcoffee_align(
        await fixture(SEQUENCES), prefix="test", compress=True, args="-out_lib=sample_lib1.tc_lib"
    )
    await check_contains(result.alignment, "test.aln.gz")
    assert result.lib is not None and result.lib.name == "sample_lib1.tc_lib", f"library: {result.lib}"


# alncompare

ALNCOMPARE_MD5 = "c77aceec520beb56f08c342e01c56a14"


@env.task
async def test_tcoffee_alncompare_seatoxin() -> None:
    scores = await tcoffee_alncompare(await fixture(SEQUENCES), await fixture(REFERENCE_MSA), meta={"id": "test"})
    assert scores.name == "setoxin-ref.scores", f"output named {scores.name}"
    await assert_md5(scores, ALNCOMPARE_MD5, label="setoxin-ref.scores")


@env.task
async def test_tcoffee_alncompare_seatoxin_compressed() -> None:
    msa = await gzipped(await fixture(SEQUENCES))
    scores = await tcoffee_alncompare(msa, await fixture(REFERENCE_MSA), meta={"id": "test"})
    assert scores.name == "setoxin-ref.fa.scores", f"output named {scores.name}"
    await assert_md5(scores, ALNCOMPARE_MD5, label="setoxin-ref.fa.scores")


# consensus


async def consensus_inputs() -> list[File]:
    sequences = await fixture(SEQUENCES)
    return [await famsa_align(sequences, prefix="test"), await magus_align(sequences, prefix="test2")]


@env.task
async def test_tcoffee_consensus_no_tree_uncompressed() -> None:
    result = await tcoffee_consensus(await consensus_inputs(), prefix="consensus", args="-output fasta_aln")
    await check_contains(result.alignment, "consensus.aln")


@env.task
async def test_tcoffee_consensus_evaluation_no_tree_uncompressed() -> None:
    result = await tcoffee_consensus(
        await consensus_inputs(), prefix="consensus", args="-output fasta_aln, score_html, sp_ascii"
    )
    await check_contains(result.alignment, "consensus.aln")


@env.task
async def test_tcoffee_consensus_tree_compressed() -> None:
    tree = await clustalo_guidetree(await fixture(SEQUENCES), prefix="test")
    result = await tcoffee_consensus(
        await consensus_inputs(), tree=tree, prefix="consensus", compress=True, args="-output fasta_aln"
    )
    await check_contains(result.alignment, "consensus.aln.gz")


# extractfrompdb


@env.task
async def test_tcoffee_extractfrompdb_seatoxin() -> None:
    pdb = (await seatoxin_structures())[0]
    out = await tcoffee_extractfrompdb(pdb, prefix="test")
    assert out.name == "test.pdb", f"output named {out.name}"
    await assert_md5(out, "d71e4feec5e11e98f9e3beb1f3ba0085", label="test.pdb")


# irmsd


@env.task
async def test_tcoffee_irmsd_seatoxin() -> None:
    irmsd = await tcoffee_irmsd(await fixture(REFERENCE_MSA), await fixture(TEMPLATE), await seatoxin_structures())
    await check_contains(irmsd, "setoxin.irmsd")


@env.task
async def test_tcoffee_irmsd_seatoxin_compressed() -> None:
    msa = await gzipped(await fixture(REFERENCE_MSA))
    irmsd = await tcoffee_irmsd(msa, await fixture(TEMPLATE), await seatoxin_structures())
    await check_contains(irmsd, "setoxin.ref.irmsd")


# regressive


@env.task
async def test_tcoffee_regressive_sequence_compressed() -> None:
    out = await tcoffee_regressive(await fixture(SEQUENCES), prefix="test", compress=True)
    await check_contains(out, "test.aln.gz")


@env.task
async def test_tcoffee_regressive_sequence_uncompressed() -> None:
    out = await tcoffee_regressive(await fixture(SEQUENCES), prefix="test")
    await check_contains(out, "test.aln")


@env.task
async def test_tcoffee_regressive_with_guide_tree() -> None:
    sequences = await fixture(SEQUENCES)
    tree = await famsa_guidetree(sequences, prefix="test")
    out = await tcoffee_regressive(sequences, tree=tree, prefix="test", compress=True)
    await check_contains(out, "test.aln.gz")


# seqreformat


@env.task
async def test_tcoffee_seqreformat_seatoxin() -> None:
    # Upstream's test doesn't load its config (-output=sim_idscore), so its
    # snapshot is of seq_reformat's default output, the sequences as FASTA.
    out = await tcoffee_seqreformat(await fixture(SEQUENCES), prefix="test")
    assert out.name == "test.txt", f"output named {out.name}"
    await assert_md5(out, "fcd4691daf120c88ec5def7ac06fb562", label="test.txt")


# tcs


async def check_tcs(result) -> None:
    await check_contains(result.tcs, "test.tcs")
    assert result.scores.name == "test.scores", f"scores named {result.scores.name}"
    lines = (await text_of(result.scores)).splitlines()
    assert lines[0] == "id,TCS" and lines[1].startswith("test,") and len(lines[1]) > len("test,"), f"scores: {lines}"


@env.task
async def test_tcoffee_tcs_seatoxin() -> None:
    await check_tcs(await tcoffee_tcs(await fixture(REFERENCE_MSA), prefix="test"))


@env.task
async def test_tcoffee_tcs_seatoxin_compressed() -> None:
    msa = await gzipped(await fixture(REFERENCE_MSA))
    await check_tcs(await tcoffee_tcs(msa, prefix="test"))


@env.task
async def test_tcoffee_tcs_seatoxin_lib() -> None:
    aligned = await tcoffee_align(
        await fixture(SEQUENCES), prefix="test", compress=True, args="-out_lib sample_lib1.tc_lib"
    )
    assert aligned.lib is not None, "tcoffee_align wrote no library"
    result = await tcoffee_tcs(await fixture(REFERENCE_MSA), lib=aligned.lib, prefix="test")
    await check_tcs(result)
    await assert_nonempty(result.tcs, label="test.tcs")


tests = [
    test_tcoffee_align_sequence_compressed,
    test_tcoffee_align_sequence_uncompressed,
    test_tcoffee_align_with_guide_tree,
    test_tcoffee_align_with_structure,
    test_tcoffee_align_with_lib,
    test_tcoffee_alncompare_seatoxin,
    test_tcoffee_alncompare_seatoxin_compressed,
    test_tcoffee_consensus_no_tree_uncompressed,
    test_tcoffee_consensus_evaluation_no_tree_uncompressed,
    test_tcoffee_consensus_tree_compressed,
    test_tcoffee_extractfrompdb_seatoxin,
    test_tcoffee_irmsd_seatoxin,
    test_tcoffee_irmsd_seatoxin_compressed,
    test_tcoffee_regressive_sequence_compressed,
    test_tcoffee_regressive_sequence_uncompressed,
    test_tcoffee_regressive_with_guide_tree,
    test_tcoffee_seqreformat_seatoxin,
    test_tcoffee_tcs_seatoxin,
    test_tcoffee_tcs_seatoxin_compressed,
    test_tcoffee_tcs_seatoxin_lib,
]
