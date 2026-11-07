"""Tests for flyte_bio.modules.trinity.

The cases and expected values are upstream's trinity tests (modules commit
9a1cc9f, Trinity 2.15.2). Upstream snapshots only the version, and checks
that the log reports success and that the assembly holds one given
transcript header; so do these.

As upstream's test config, Trinity runs on 1 CPU (assemblies can change with
the thread count), with upstream's 15 GB of test memory, so ``--max_memory``
is 12G as there.

Upstream's "multiple" cases pass the same files twice. Trinity's reads are
staged side by side under their names, so here the repeats are renamed copies
of the same files (Trinity reads their content, not their names).
"""

import gzip
import shutil
import tempfile
from pathlib import Path

import flyte
from flyte.io import File

from flyte_bio.modules.trinity import TrinityResult, read_trinity_results, trinity, trinity_cmd
from tests.framework import env, fixture

FASTQ = "genomics/homo_sapiens/illumina/fastq/"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def gunzipped_text(file: File) -> str:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode()


async def renamed(file: File, name: str) -> File:
    """A copy of ``file`` named ``name``."""
    local = await file.download()
    dest = Path(tempfile.mkdtemp()) / name
    shutil.copy(local, dest)
    return await File.from_local(str(dest))


async def run_trinity(reads: list[File], single_end: bool) -> TrinityResult:
    """Run Trinity as upstream's test does: 1 CPU, 15 GB."""
    names = [r.name for r in reads]
    if single_end:
        left, right = names, []
    else:
        left, right = names[0::2], names[1::2]
    # The cache key ignores resources, so salt it to keep this run apart from a default one.
    task = trinity_cmd.as_task().override(
        resources=flyte.Resources(cpu=1, memory="15Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=1,mem=15Gi")
    )
    results = await task(reads=reads, left=",".join(left), right=",".join(right), seq_type="fq", prefix="test", args="")
    return await read_trinity_results(results, "test")


async def check(result: TrinityResult, header: str) -> None:
    assert "All commands completed successfully." in await text(result.log), "trinity log lacks success"
    fasta = await gunzipped_text(result.transcript_fasta)
    assert header in fasta.splitlines(), f"assembly lacks {header!r}"


@env.task
async def test_trinity_paired_end() -> None:
    reads = [await fixture(FASTQ + "test_rnaseq_1.fastq.gz"), await fixture(FASTQ + "test_rnaseq_2.fastq.gz")]
    await check(await run_trinity(reads, single_end=False), ">TRINITY_DN7_c0_g1_i1 len=299 path=[0:0-151 2:152-298]")


@env.task
async def test_trinity_single_end() -> None:
    reads = [await fixture(FASTQ + "test_rnaseq_1.fastq.gz")]
    await check(await run_trinity(reads, single_end=True), ">TRINITY_DN10_c0_g1_i1 len=534 path=[0:0-533]")


@env.task
async def test_trinity_paired_end_multiple() -> None:
    r1 = await fixture(FASTQ + "test_rnaseq_1.fastq.gz")
    r2 = await fixture(FASTQ + "test_rnaseq_2.fastq.gz")
    reads = [r1, r2, await renamed(r1, "copy_test_rnaseq_1.fastq.gz"), await renamed(r2, "copy_test_rnaseq_2.fastq.gz")]
    await check(await run_trinity(reads, single_end=False), ">TRINITY_DN12_c0_g1_i1 len=225 path=[0:0-224]")


@env.task
async def test_trinity_single_end_multiple() -> None:
    r1 = await fixture(FASTQ + "test_rnaseq_1.fastq.gz")
    reads = [r1, await renamed(r1, "copy_test_rnaseq_1.fastq.gz")]
    await check(await run_trinity(reads, single_end=True), ">TRINITY_DN12_c0_g1_i1 len=268 path=[0:0-267]")


@env.task
async def test_trinity_wrapper() -> None:
    """The wrapper's argument handling: rejects odd paired input and repeated names."""
    r1 = await fixture(FASTQ + "test_rnaseq_1.fastq.gz")
    for reads, single_end in [([r1], False), ([r1, r1], True)]:
        try:
            await trinity(reads, single_end=single_end, prefix="test")
        except ValueError:
            continue
        raise AssertionError(f"trinity accepted {[r.name for r in reads]} (single_end={single_end})")


tests = [
    test_trinity_paired_end,
    test_trinity_single_end,
    test_trinity_paired_end_multiple,
    test_trinity_single_end_multiple,
    test_trinity_wrapper,
]
