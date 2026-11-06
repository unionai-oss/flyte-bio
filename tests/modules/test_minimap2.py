"""Tests for flyte_bio.modules.minimap2.

The cases and expected values are upstream's minimap2/index and
minimap2/align tests (modules commit 9a1cc9f, minimap2 2.30-r1287, samtools
1.23.1). minimap2/index's snapshot is the ``.mmi``'s md5. minimap2/align's
snapshots are the BAM header and the md5 of the read sequences (nft-bam
``getReadsMD5``). The header's ``@PG`` lines record the thread counts, which
upstream's tests ran with 2 CPUs (``-t 2``, ``-@ 1``); the counts are
normalised to those before comparing, since they follow the task's CPUs.
"""

import asyncio
import hashlib
import re

from flyte.io import File

from flyte_bio.modules.minimap2 import Minimap2AlignResult, minimap2_align, minimap2_index
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_md5, assert_reads_md5, env, fixture

FASTA = "genomics/sarscov2/genome/genome.fasta"
FASTQ = "genomics/sarscov2/illumina/fastq/"
BAM = "genomics/homo_sapiens/illumina/bam/test3.single_end.markduplicates.sorted.bam"

HD = "@HD\tVN:1.6\tSO:coordinate"
SARSCOV2_SQ = ["@SQ\tSN:MT192765.1\tLN:29829"]
MINIMAP2_PG = "@PG\tID:minimap2\tPN:minimap2\tVN:2.30-r1287\tCL:minimap2 -t 2 -a "
SORT_PG = "@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.23.1\tCL:samtools sort -@ 1 -o test.bam"
SORT_INDEX_PG = SORT_PG + "##idx##test.bam.bai --write-index"


async def header(bam: File) -> list[str]:
    """``bam``'s header lines, with the @PG thread counts set to upstream's."""
    text = await samtools_view(bam=bam, args="-H --no-PG")
    async with text.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    lines = [re.sub(r"CL:minimap2 -t \d+ ", "CL:minimap2 -t 2 ", line) for line in lines]
    return [re.sub(r"CL:samtools sort -@ \d+ ", "CL:samtools sort -@ 1 ", line) for line in lines]


async def check_bam(
    result: Minimap2AlignResult, expected_header: list[str], reads_md5: str, label: str, index: bool = False
) -> None:
    assert result.bam is not None and result.paf is None, result
    assert result.bam.path.endswith("/test.bam"), result.bam.path
    actual = await header(result.bam)
    assert actual == expected_header, f"{label}: header {actual} != {expected_header}"
    await assert_reads_md5(result.bam, reads_md5, label=label)
    if index:
        assert result.index is not None and result.index.path.endswith("/test.bam.bai"), result.index
    else:
        assert result.index is None, result.index


@env.task
async def test_minimap2_index() -> None:
    # upstream case: minimap2/index "minimap2 index"
    index = await minimap2_index(await fixture(FASTA))
    assert index.path.endswith("/genome.mmi"), index.path
    await assert_md5(index, "72e450f12dc691e763c697463bdb1571", label="genome.mmi")


@env.task
async def test_minimap2_align_fastq() -> None:
    # upstream cases: "sarscov2 - fastq, fasta, true, [], false, false" and "..., true, 'bai', false, false"
    fasta, reads = await asyncio.gather(fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"))
    plain, indexed = await asyncio.gather(
        minimap2_align([reads], fasta, prefix="test"),
        minimap2_align([reads], fasta, prefix="test", bam_index_extension="bai"),
    )
    pg = MINIMAP2_PG + "genome.fasta test_1.fastq.gz"
    md5 = "f194745c0ccfcb2a9c0aee094a08750"
    await check_bam(plain, [HD, *SARSCOV2_SQ, pg, SORT_PG], md5, "fastq")
    await check_bam(indexed, [HD, *SARSCOV2_SQ, pg, SORT_INDEX_PG], md5, "fastq, bai", index=True)


@env.task
async def test_minimap2_align_paired_fastq() -> None:
    # upstream case: "sarscov2 - [fastq1, fastq2], fasta, true, false, false"
    fasta, r1, r2 = await asyncio.gather(
        fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    result = await minimap2_align([r1, r2], fasta, prefix="test")
    pg = MINIMAP2_PG + "genome.fasta test_1.fastq.gz test_2.fastq.gz"
    await check_bam(result, [HD, *SARSCOV2_SQ, pg, SORT_PG], "1bc392244f228bf52cf0b5a8f6a654c9", "paired fastq")


@env.task
async def test_minimap2_align_fastq_no_reference() -> None:
    # upstream case: "sarscov2 - fastq, [], true, false, false": the reads are aligned to themselves,
    # so the header has one @SQ per read (100), checked by their count and md5.
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await minimap2_align([reads], None, prefix="test")
    assert result.bam is not None, result
    lines = await header(result.bam)
    sq = [line for line in lines if line.startswith("@SQ")]
    assert len(sq) == 100, f"{len(sq)} @SQ lines, expected 100"
    assert sq[0] == "@SQ\tSN:ERR5069949.2151832\tLN:150" and sq[-1] == "@SQ\tSN:ERR5069949.2415814\tLN:150", sq
    sq_md5 = hashlib.md5("\n".join(sq).encode()).hexdigest()
    assert sq_md5 == "3696d105433855a04c6aff5760cbf968", f"@SQ lines md5 {sq_md5}"
    expected = [HD, *sq, MINIMAP2_PG + "test_1.fastq.gz test_1.fastq.gz", SORT_PG]
    assert lines == expected, f"header {lines[:1] + lines[-2:]} != {expected[:1] + expected[-2:]}"
    await assert_reads_md5(result.bam, "16c1c651f8ec67383bcdee3c55aed94f", label="fastq, no reference")


@env.task
async def test_minimap2_align_bam() -> None:
    # upstream cases: "sarscov2 - bam, fasta, true, [], false, false" and "..., true, 'bai', false, false"
    fasta, bam = await asyncio.gather(fixture(FASTA), fixture(BAM))
    plain, indexed = await asyncio.gather(
        minimap2_align([bam], fasta, prefix="test"),
        minimap2_align([bam], fasta, prefix="test", bam_index_extension="bai"),
    )
    pg = MINIMAP2_PG + "genome.fasta -"
    md5 = "5d426b9a5f5b2c54f1d7f1e4c238ae94"
    await check_bam(plain, [HD, *SARSCOV2_SQ, pg, SORT_PG], md5, "bam")
    await check_bam(indexed, [HD, *SARSCOV2_SQ, pg, SORT_INDEX_PG], md5, "bam, bai", index=True)


@env.task
async def test_minimap2_align_bam_no_reference() -> None:
    # upstream case: "sarscov2 - bam, [], true, false, false" must fail: alignment input needs a reference.
    bam = await fixture(BAM)
    try:
        await minimap2_align([bam], None, prefix="test")
    except ValueError as e:
        assert "needs a reference" in str(e), e
    else:
        raise AssertionError("minimap2_align accepted alignment input without a reference")


@env.task
async def test_minimap2_align_paf_and_index() -> None:
    # Not an upstream case (upstream tests PAF output and .mmi references only as stubs): align to
    # the built index, which must give the same reads as the FASTA, and write PAF with CIGARs.
    fasta, reads = await asyncio.gather(fixture(FASTA), fixture(FASTQ + "test_1.fastq.gz"))
    index = await minimap2_index(fasta)
    bam, paf = await asyncio.gather(
        minimap2_align([reads], index, prefix="test"),
        minimap2_align([reads], fasta, prefix="test", bam_format=False, cigar_paf_format=True),
    )
    assert bam.bam is not None, bam
    await assert_reads_md5(bam.bam, "f194745c0ccfcb2a9c0aee094a08750", label="fastq, mmi")
    assert paf.bam is None and paf.index is None and paf.paf is not None, paf
    assert paf.paf.path.endswith("/test.paf"), paf.paf.path
    async with paf.paf.open("rb") as fh:
        records = [line.split("\t") for line in bytes(await fh.read()).decode().splitlines()]
    assert records, "empty PAF"
    for fields in records:
        assert len(fields) >= 12 and fields[5] == "MT192765.1", fields[:12]
        assert any(f.startswith("cg:Z:") for f in fields[12:]), f"no CIGAR in {fields[:12]}"


tests = [
    test_minimap2_index,
    test_minimap2_align_fastq,
    test_minimap2_align_paired_fastq,
    test_minimap2_align_fastq_no_reference,
    test_minimap2_align_bam,
    test_minimap2_align_bam_no_reference,
    test_minimap2_align_paf_and_index,
]
