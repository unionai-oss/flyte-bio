"""Tests for flyte_bio.modules.yara.

The cases and expected values are upstream's yara/index and yara/mapper
tests (modules commit 9a1cc9f, Yara 1.0.2, samtools 1.16.1). yara/index
snapshots the index files' md5s. yara/mapper snapshots each BAM's name,
read md5 (``getReadsMD5``) and header md5. The header records Yara's and
samtools' command lines, so it's checked by content instead: the reference,
Yara's program record and the coordinate sort. As upstream, the mapper runs
on one CPU, since Yara's read order depends on its thread count.
"""

import asyncio

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.samtools import samtools_view
from flyte_bio.modules.yara import YaraMapperResult, yara_index, yara_mapper, yara_mapper_cmd
from tests.framework import assert_md5, assert_reads_md5, env, fixture

FASTA = "genomics/sarscov2/genome/genome.fasta"
FASTQ = "genomics/sarscov2/illumina/fastq/"
INDEX_MD5S = {
    "genome.fasta.lf.drp": "3ef99a87a4e44513f46d42f4261f7842",
    "genome.fasta.lf.drs": "55a54008ad1ba589aa210d2629c1df41",
    "genome.fasta.lf.drv": "cf6408307fe9fd7f99c33f521bf95550",
    "genome.fasta.lf.pst": "e8daba34298e99e42942435286f9b3f0",
    "genome.fasta.rid.concat": "1e4e4c88ddeaf907a12f02f0d88367c5",
    "genome.fasta.rid.limits": "8b814661f30a0c9e350bfbcb454930ce",
    "genome.fasta.sa.ind": "464314583efb5f07260b0efecc29a1ce",
    "genome.fasta.sa.len": "45677f66c28c79c02250ceb8b58645e8",
    "genome.fasta.sa.val": "ce57cc82e2d3ae7b9824210f54168ce9",
    "genome.fasta.txt.concat": "6074d1933c9e7e5ab05fa0def5ce28c0",
    "genome.fasta.txt.limits": "4480a068db603e4c9a27bc4fa9ceaf14",
    "genome.fasta.txt.size": "063987b3c3f747be7d2b8043c9d91000",
}


async def header_lines(bam: File) -> list[str]:
    header = await samtools_view(bam=bam, args="-H")
    async with header.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def check_bams(result: YaraMapperResult, expected: dict[str, str], label: str) -> None:
    names = [bam.name for bam in result.bams]
    assert names == list(expected), f"{label}: BAMs are {names}, expected {list(expected)}"
    assert [bai.name for bai in result.bais] == [f"{name}.bai" for name in names], f"{label}: indexes"
    for bam, md5 in zip(result.bams, expected.values(), strict=True):
        await assert_reads_md5(bam, md5, label=f"{label} {bam.name}")
        header = await header_lines(bam)
        assert "@SQ\tSN:MT192765.1\tLN:29829" in header, f"{label} {bam.name}: {header}"
        assert any(line.startswith("@PG") and "yara_mapper" in line for line in header), f"{label}: {header}"
        assert any(line.startswith("@HD") and "SO:coordinate" in line for line in header), f"{label}: {header}"


async def map_on_one_cpu(reads_1: File, reads_2: File | None, index: Dir) -> YaraMapperResult:
    # The cache key ignores resources, so salt it to keep the 1-CPU run apart from a default one.
    task = yara_mapper_cmd.override(
        resources=flyte.Resources(cpu=1, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=1")
    )
    results = await task(reads_1=reads_1, reads_2=reads_2, index=index, prefix="test", args="")
    names = ["test.mapped.bam"] if reads_2 is None else ["test_1.mapped.bam", "test_2.mapped.bam"]
    bams = [await results.get_file(name) for name in names]
    bais = [await results.get_file(f"{name}.bai") for name in names]
    assert all(bams) and all(bais), f"yara mapper wrote {[f.name for f in await results.list_files()]}"
    return YaraMapperResult(bams=bams, bais=bais, results=results)


async def sarscov2_index() -> Dir:
    return await yara_index(await fixture(FASTA))


@env.task
async def test_yara_index() -> None:
    # upstream case: yara/index "sarscov2 - bam"
    index = await sarscov2_index()
    files = sorted(f.name for f in await index.list_files())
    assert files == sorted(INDEX_MD5S), f"index files {files}"
    for name, md5 in INDEX_MD5S.items():
        await assert_md5(await index.get_file(name), md5, label=name)


@env.task
async def test_yara_mapper_single_end() -> None:
    # upstream case: yara/mapper "sarscov2 - single_end"
    index, reads = await asyncio.gather(sarscov2_index(), fixture(FASTQ + "test_1.fastq.gz"))
    result = await map_on_one_cpu(reads, None, index)
    await check_bams(result, {"test.mapped.bam": "c0753fc123f2bce542a2e425b397c478"}, "single-end")


@env.task
async def test_yara_mapper_paired_end() -> None:
    # upstream case: yara/mapper "sarscov2 - paired_end"
    index, r1, r2 = await asyncio.gather(
        sarscov2_index(), fixture(FASTQ + "test_1.fastq.gz"), fixture(FASTQ + "test_2.fastq.gz")
    )
    result = await map_on_one_cpu(r1, r2, index)
    expected = {
        "test_1.mapped.bam": "dc3035ba9a7864563016636d5145726b",
        "test_2.mapped.bam": "4d493e899be3833d7b3d8c111c0a740c",
    }
    await check_bams(result, expected, "paired-end")


@env.task
async def test_yara_mapper_wrapper() -> None:
    # Not an upstream case: the wrapper at the default CPUs finds every output.
    index, reads = await asyncio.gather(sarscov2_index(), fixture(FASTQ + "test_1.fastq.gz"))
    result = await yara_mapper(reads, None, index, "test")
    assert [bam.name for bam in result.bams] == ["test.mapped.bam"], f"BAMs are {[b.name for b in result.bams]}"


tests = [test_yara_index, test_yara_mapper_single_end, test_yara_mapper_paired_end, test_yara_mapper_wrapper]
