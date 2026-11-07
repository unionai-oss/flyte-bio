"""Tests for flyte_bio.modules.modkit.

The cases and expected values are upstream's modkit/callmods, modkit/pileup
and modkit/repair tests (modules commit 9a1cc9f, modkit 0.6.4). call-mods and
repair snapshot the md5 of the output BAM's read sequences (nft-bam
``getReadsMD5``) and the log's name; pileup snapshots the bedMethyl files'
md5s, which nf-test computes decompressed, and no log.

Upstream's repair test runs on 2 CPUs, since modkit repair writes reads out
of order with 4 or more threads; so does this one (a cache salt keeps the
2-CPU run apart from a default one). The pileup BED is upstream's
one-line ``chr22.bed``.
"""

import asyncio
import tempfile
from pathlib import Path

import flyte
from flyte.io import File

from flyte_bio.modules.modkit import (
    ModkitBamResult,
    ModkitPileupResult,
    bam_and_log,
    modkit_callmods,
    modkit_pileup,
    modkit_repair_cmd,
)
from tests.framework import assert_gunzipped_md5, assert_reads_md5, env, fixture

NANOPORE = "genomics/homo_sapiens/nanopore/bam/"
PHASED_BAM = NANOPORE + "test.sorted.phased.bam"
GENOME = "genomics/homo_sapiens/genome/genome.fasta"
ARABIDOPSIS = "genomics/eukaryotes/arabidopsis_thaliana/genome/"
MODIFIED_BASES = "--modified-bases 5mC 5hmC"


async def check_bam_and_log(result: ModkitBamResult, reads_md5: str, label: str) -> None:
    assert result.bam.name == "test.bam", result.bam.name
    assert result.log.name == "test.log", result.log.name
    await assert_reads_md5(result.bam, reads_md5, label=label)


async def check_pileup(result: ModkitPileupResult, expected: dict[str, str]) -> None:
    names = [f.name for f in result.beds]
    assert names == sorted(expected), f"bedMethyl files {names} != {sorted(expected)}"
    for bed in result.beds:
        await assert_gunzipped_md5(bed, expected[bed.name], label=bed.name)
    assert result.log is None, f"unexpected log {result.log}"


async def phased_inputs() -> tuple[File, File]:
    return await asyncio.gather(fixture(PHASED_BAM), fixture(PHASED_BAM + ".bai"))


async def reference() -> tuple[File, File]:
    return await asyncio.gather(fixture(GENOME), fixture(GENOME + ".fai"))


async def chr22_bed() -> File:
    path = Path(tempfile.mkdtemp()) / "chr22.bed"
    path.write_text("chr22\t0\t1000\n")
    return await File.from_local(str(path))


@env.task
async def test_modkit_callmods() -> None:
    # upstream case: modkit/callmods "homo_sapiens - bam"
    bam = await fixture(NANOPORE + "bc_anchored_10_reads.sorted.bam")
    result = await modkit_callmods(bam, prefix="test")
    await check_bam_and_log(result, "4e7b32d48deea26eae8893d0de451f73", "call-mods")


@env.task
async def test_modkit_pileup_no_reference() -> None:
    # upstream case: modkit/pileup "[bam, bai], [], []"
    bam, bai = await phased_inputs()
    result = await modkit_pileup(bam, bai, prefix="test")
    await check_pileup(result, {"test.bed.gz": "1ec207494e615704e7806d4a6621de93"})


@env.task
async def test_modkit_pileup_reference() -> None:
    # upstream case: modkit/pileup "[bam, bai], [fasta, fai], []"
    (bam, bai), (fasta, fai) = await asyncio.gather(phased_inputs(), reference())
    result = await modkit_pileup(bam, bai, fasta, fai, prefix="test")
    await check_pileup(result, {"test.bed.gz": "1ec207494e615704e7806d4a6621de93"})


@env.task
async def test_modkit_pileup_reference_bed() -> None:
    # upstream case: modkit/pileup "[bam, bai], [fasta, fai], bed", with --modified-bases 5mC 5hmC
    (bam, bai), (fasta, fai), bed = await asyncio.gather(phased_inputs(), reference(), chr22_bed())
    result = await modkit_pileup(bam, bai, fasta, fai, bed, prefix="test", args=MODIFIED_BASES)
    await check_pileup(result, {"test.bed.gz": "cd173b78a199424358268f7f15fc1878"})


@env.task
async def test_modkit_pileup_phased() -> None:
    # upstream case: modkit/pileup "[bam, bai], [fasta, fai], [], phased"
    (bam, bai), (fasta, fai) = await asyncio.gather(phased_inputs(), reference())
    result = await modkit_pileup(bam, bai, fasta, fai, prefix="test", args="--phased " + MODIFIED_BASES)
    await check_pileup(
        result,
        {
            "test_combined.bed.gz": "32452a02e350e5ac558ab68bbbaccc85",
            "test_hp1.bed.gz": "abd206b393ead3f3a574a4c34dc7429a",
            "test_hp2.bed.gz": "5092dde96578433a30c4a7c86f3cd94d",
        },
    )


@env.task
async def test_modkit_pileup_phased_bed() -> None:
    # upstream case: modkit/pileup "[bam, bai], [fasta, fai], bed, phased"
    (bam, bai), (fasta, fai), bed = await asyncio.gather(phased_inputs(), reference(), chr22_bed())
    result = await modkit_pileup(bam, bai, fasta, fai, bed, prefix="test", args="--phased " + MODIFIED_BASES)
    await check_pileup(
        result,
        {
            "test_combined.bed.gz": "cd173b78a199424358268f7f15fc1878",
            "test_hp1.bed.gz": "d67280b345e95029cb5c0ea771128ba1",
            "test_hp2.bed.gz": "6438ea42b16a21db197ddf73ec665685",
        },
    )


@env.task
async def test_modkit_repair() -> None:
    # upstream case: modkit/repair "test - donor_bam - acceptor_bam", on 2 CPUs as upstream's test config
    donor, acceptor = await asyncio.gather(
        fixture(ARABIDOPSIS + "ont_Col-0_tets_data_Chr1_5mb_7mb.bam"),
        fixture(ARABIDOPSIS + "trimmed_ont_Col-0_tets_data_Chr1_5mb_7mb.bam"),
    )
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    task = modkit_repair_cmd.as_task().override(
        resources=flyte.Resources(cpu=2, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=2")
    )
    results = await task(donor_bam=donor, acceptor_bam=acceptor, prefix="test", args="")
    result = await bam_and_log(results, "test", "repair")
    await check_bam_and_log(result, "4a25b01c940d7ee25fc8b55763c4b899", "repair")


tests = [
    test_modkit_callmods,
    test_modkit_pileup_no_reference,
    test_modkit_pileup_reference,
    test_modkit_pileup_reference_bed,
    test_modkit_pileup_phased,
    test_modkit_pileup_phased_bed,
    test_modkit_repair,
]
