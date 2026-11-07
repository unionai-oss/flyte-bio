"""Tests for flyte_bio.modules.sortmerna.

The cases are upstream's sortmerna tests (modules commit 9a1cc9f, SortMeRNA
4.3.7). Upstream snapshots the index files' md5s (but not the ``.stats``
file, which records paths and times), the names of the filtered reads and
log, and a line of the log.
"""

from flyte.io import Dir, File

from flyte_bio.modules.sortmerna import SortmernaResult, sortmerna
from tests.framework import assert_md5, assert_nonempty, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
INDEX_MD5S = {
    "2415186086593376314.bursttrie_0.dat": "74f7f020e8d46e24a8a2e9c5fbcd564a",
    "2415186086593376314.kmer_0.dat": "4a0bcb71b120f6a6949b7969292ef2e7",
    "2415186086593376314.pos_0.dat": "bc2875e4cc4017707306565e396839ef",
}


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def references() -> list[File]:
    return [await fixture("genomics/sarscov2/genome/genome.fasta")]


async def check_index(index: Dir, label: str) -> None:
    files = {f.name: f for f in await index.list_files() if not f.name.endswith(".stats")}
    assert sorted(files) == sorted(INDEX_MD5S), f"{label}: index holds {sorted(files)}"
    for name, md5 in INDEX_MD5S.items():
        await assert_md5(files[name], md5, label=f"{label} {name}")


async def check_filtered(result: SortmernaResult, names: list[str], passing: str, label: str) -> None:
    assert [f.name for f in result.reads] == names, f"{label}: reads {[f.name for f in result.reads]}"
    for f in result.reads:
        await assert_nonempty(f, label=f"{label} {f.name}")
    assert result.log is not None and result.log.name == "test.sortmerna.log", f"{label}: log {result.log}"
    line = f"Total reads passing E-value threshold = {passing} (100.00)"
    assert line in await text(result.log), f"{label}: log lacks {line!r}"
    await check_index(result.index, label)


@env.task
async def test_sortmerna_indexing_only() -> None:
    # upstream case: sortmerna "sarscov2 indexing only" (ext.args --index 1)
    result = await sortmerna(await references(), prefix="test", args="--index 1")
    assert result.reads == [] and result.log is None
    await check_index(result.index, "indexing only")


@env.task
async def test_sortmerna_single_end() -> None:
    # upstream case: sortmerna "sarscov2 single_end"
    result = await sortmerna(await references(), await fixture(FASTQ + "test_1.fastq.gz"), prefix="test")
    await check_filtered(result, ["test.non_rRNA.fastq.gz"], "100", "single-end")


@env.task
async def test_sortmerna_paired_end() -> None:
    # upstream case: sortmerna "sarscov2 paired_end"
    r1, r2 = await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")
    result = await sortmerna(await references(), r1, r2, prefix="test")
    await check_filtered(result, ["test_1.non_rRNA.fastq.gz", "test_2.non_rRNA.fastq.gz"], "200", "paired-end")


@env.task
async def test_sortmerna_single_end_premade_index() -> None:
    # upstream case: sortmerna "sarscov2 single_end premade_index" (index --index 1, then --index 0)
    refs = await references()
    built = await sortmerna(refs, prefix="test2", args="--index 1")
    reads = await fixture(FASTQ + "test_1.fastq.gz")
    result = await sortmerna(refs, reads, index=built.index, prefix="test", args="--index 0")
    await check_filtered(result, ["test.non_rRNA.fastq.gz"], "100", "premade index")


tests = [
    test_sortmerna_indexing_only,
    test_sortmerna_single_end,
    test_sortmerna_paired_end,
    test_sortmerna_single_end_premade_index,
]
