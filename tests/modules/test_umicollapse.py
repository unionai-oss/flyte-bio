"""Tests for flyte_bio.modules.umicollapse.

The cases are upstream's umicollapse tests (modules commit 9a1cc9f,
UMICollapse 1.1.0; ext.prefix ``<id>.dedup``). The FASTQ cases assert
upstream's md5s. For BAMs upstream snapshots ``getSamLinesMD5``, which hashes
htsjdk's SAM text including the ``@PG`` command line (with the jar's path), so
it can't be reproduced; those cases check the header and invariants of
deduplication instead: the output's reads are a subset of the input's (fewer,
for the single-end BAM), and every mapped position of the input keeps at least
one read.

Upstream's paired-end FASTQ case stages both mates but runs UMICollapse on
mate 1 only (its command reads one input), so the port does the same.
"""

from flyte.io import File

from flyte_bio.modules.samtools import samtools_view
from flyte_bio.modules.umicollapse import umicollapse
from tests.framework import assert_gunzipped_md5, assert_nonempty, env, fixture

BAM = "genomics/sarscov2/illumina/bam/"
FASTQ = "genomics/sarscov2/illumina/fastq/"


async def sam_lines(bam: File) -> list[str]:
    sam = await samtools_view(bam=bam, args="-h")
    async with sam.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def check_dedup(name: str, removes_reads: bool) -> None:
    bam = await fixture(BAM + f"{name}.bam")
    bai = await fixture(BAM + f"{name}.bam.bai")
    result = await umicollapse(bam, "bam", bai=bai, prefix="test.dedup")
    assert result.bam is not None and result.bam.name == "test.dedup.bam", result.bam
    assert result.log.name == "test.dedup_UMICollapse.log", result.log.name
    await assert_nonempty(result.log, label=f"{name} log")

    before = [line.split("\t") for line in await sam_lines(bam) if not line.startswith("@")]
    after_lines = await sam_lines(result.bam)
    assert any(line.startswith("@HD") for line in after_lines), f"{name}: output BAM lacks its header"
    after = [line.split("\t") for line in after_lines if not line.startswith("@")]
    assert after, f"{name}: no reads left after deduplication"
    assert len(after) <= len(before), f"{name}: {len(after)} reads left of {len(before)}"
    if removes_reads:
        assert len(after) < len(before), f"{name}: {len(after)} reads left of {len(before)}, none removed"
    records = {"\t".join(r[:11]) for r in before}
    missing = [r[0] for r in after if "\t".join(r[:11]) not in records]
    assert not missing, f"{name}: reads not in the input: {missing[:5]}"

    def positions(rows: list[list[str]]) -> set[tuple[str, str, str]]:
        # (reference, position, reverse strand) of mapped reads
        return {(r[2], r[3], str(int(r[1]) & 16)) for r in rows if not int(r[1]) & 4}

    lost = positions(before) - positions(after)
    assert not lost, f"{name}: positions left without a read: {sorted(lost)[:5]}"


@env.task
async def test_umicollapse_bam_single_end() -> None:
    # upstream case: umicollapse "umicollapse single end test"
    await check_dedup("test.single_end.umi.sorted", removes_reads=True)


@env.task
async def test_umicollapse_bam_paired_end() -> None:
    # upstream case: umicollapse "umicollapse paired tests". Upstream passes no
    # --paired, so mates are deduplicated one by one; on this BAM none is a duplicate.
    await check_dedup("test.paired_end.umi.sorted", removes_reads=False)


@env.task
async def test_umicollapse_fastq_single_end() -> None:
    # upstream case: umicollapse "umicollapse fastq test (single-end)"
    reads = await fixture(FASTQ + "test.umi_extract_single.fastq.gz")
    result = await umicollapse(reads, "fastq", prefix="test.dedup")
    assert result.fastq is not None and result.fastq.name == "test.dedup.fastq.gz", result.fastq
    await assert_gunzipped_md5(result.fastq, "2e602ed23eb87f434e4f0a9e491c0310", label="single-end fastq")


@env.task
async def test_umicollapse_fastq_paired_end() -> None:
    # upstream case: umicollapse "umicollapse fastq test (paired-end)": deduplicates mate 1 only
    reads = await fixture(FASTQ + "test.umi_extract_1.fastq.gz")
    result = await umicollapse(reads, "fastq", prefix="test.dedup")
    assert result.fastq is not None
    await assert_gunzipped_md5(result.fastq, "721a84a2accac988d636e837c60e47bc", label="paired-end mate 1 fastq")


tests = [
    test_umicollapse_bam_single_end,
    test_umicollapse_bam_paired_end,
    test_umicollapse_fastq_single_end,
    test_umicollapse_fastq_paired_end,
]
