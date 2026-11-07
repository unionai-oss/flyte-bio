"""Tests for flyte_bio.modules.chromap.

The cases and expected values are upstream's chromap/index and
chromap/chromap tests (modules commit 9a1cc9f, chromap 0.3.2). chromap/index
snapshots only the index's name. chromap/chromap snapshots the decompressed
BED's md5 for single- and paired-end reads; its ``--SAM`` case snapshots
nothing (its ``snapshot(...)`` has no ``.match()``), so here the BAM is
checked by content: every pair's reads are there, mapped to the genome.
"""

import asyncio
from collections import Counter

from flyte_bio.modules.chromap import chromap_chromap, chromap_index
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_gunzipped_md5, env, fixture

FASTA = "genomics/sarscov2/genome/genome.fasta"
FASTQ = "genomics/sarscov2/illumina/fastq/"


@env.task
async def test_chromap_index() -> None:
    # upstream case: chromap/index "sarscov2 - fasta"
    index = await chromap_index(await fixture(FASTA))
    assert index.name == "genome.index", f"index is {index.name}, expected genome.index"


@env.task
async def test_chromap_chromap_bed() -> None:
    # upstream cases: chromap/chromap "sarscov2 - single_end" and "sarscov2 - paired_end"
    fasta, r1, r2 = await asyncio.gather(
        fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    index = await chromap_index(fasta)
    single, paired = await asyncio.gather(
        chromap_chromap([r1], fasta, index, prefix="test", single_end=True),
        chromap_chromap([r1, r2], fasta, index, prefix="test"),
    )
    for label, result, md5 in (
        ("single-end", single, "b23aa9dae812fe84c308df6e18d850be"),
        ("paired-end", paired, "f658c1e94adef751ad08e132aac92712"),
    ):
        assert result.output.name == "test.bed.gz", f"{label}: output is {result.output.name}"
        await assert_gunzipped_md5(result.output, md5, label=f"{label} test.bed.gz")


@env.task
async def test_chromap_chromap_sam() -> None:
    # upstream case: chromap/chromap "sarscov2 - paired_end - bam" (ext.args '--SAM')
    fasta, r1, r2 = await asyncio.gather(
        fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    index = await chromap_index(fasta)
    result = await chromap_chromap([r1, r2], fasta, index, prefix="test", args="--SAM")
    assert result.output.name == "test.bam", f"output is {result.output.name}, expected test.bam"
    sam = await samtools_view(bam=result.output)
    async with sam.open("rb") as fh:
        records = [line.split("\t") for line in bytes(await fh.read()).decode().splitlines()]
    assert records, "chromap wrote no alignments"
    contigs = {r[2] for r in records}
    assert contigs == {"MT192765.1"}, f"alignments on {contigs}, expected MT192765.1"
    unpaired = [name for name, count in Counter(r[0] for r in records).items() if count != 2]
    assert not unpaired, f"{len(unpaired)} reads have no mate, e.g. {unpaired[:3]}"


tests = [test_chromap_index, test_chromap_chromap_bed, test_chromap_chromap_sam]
