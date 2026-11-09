"""Tests for flyte_bio.modules.lima.

The cases and expected md5s are upstream's lima tests (modules commit
9a1cc9f, lima 2.12.0), with their test config's ``--isoseq --peek-guess`` and
prefix ``<id>.fl``. nf-test hashes every output as it is, except gzipped reads,
which it hashes decompressed. As upstream's snapshots show, a FASTA/FASTQ run
writes no ``.xml`` or ``.json`` and a BAM run no FASTA/FASTQ. The BAM
header's ``@PG`` command line records ``-j``: upstream's snapshot is of a
2-thread run, the module's default.
"""

import asyncio
import os

from flyte.io import File

from flyte_bio.modules.lima import LimaResult, lima
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

PACBIO = "genomics/homo_sapiens/pacbio/"
PRIMERS = PACBIO + "fasta/primers.fasta"
ARGS = "--isoseq --peek-guess"
PAIR = "test.fl.NEB_5p--NEB_Clontech_3p"


def names(files: list[File]) -> list[str]:
    return [os.path.basename(f.path) for f in files]


async def run(reads: str) -> LimaResult:
    ccs, primers = await asyncio.gather(fixture(PACBIO + reads), fixture(PRIMERS))
    return await lima(ccs, primers, prefix="test.fl", args=ARGS)


async def check_reports(result: LimaResult, md5s: dict[str, str]) -> None:
    for kind, md5 in md5s.items():
        found = getattr(result, kind)
        assert found is not None, f"lima wrote no test.fl.lima.{kind}"
        await assert_md5(found, md5, label=f"test.fl.lima.{kind}")


async def check_bam(result: LimaResult, bam_md5: str, pbi_md5: str) -> None:
    assert names(result.bams) == [f"{PAIR}.bam"], names(result.bams)
    assert names(result.pbis) == [f"{PAIR}.bam.pbi"], names(result.pbis)
    await assert_md5(result.bams[0], bam_md5, label=f"{PAIR}.bam")
    await assert_md5(result.pbis[0], pbi_md5, label=f"{PAIR}.bam.pbi")


# The FASTA and FASTQ inputs hold the same reads, so their outputs differ only in format.
FASTA_REPORTS = {
    "counts": "a4ceaa408be334eaa711577e95f8730e",
    "report": "bd4a8bde17471563cf91aab4c787911d",
    "summary": "03be2311ba4afb878d8e547ab38c11eb",
    "clips": "1012bc8874a14836f291bac48e8482a4",
    "guess": "651e5f2b438b8ceadb3e06a2177e1818",
}
FASTQ_REPORTS = {
    "counts": "767b687e6eda7b24cd0e577f527eb2f0",
    "report": "ad2a9b1eeb4cda4a1f69ef4b7520b5fd",
    "summary": "e91d3c386aaf4effa63f33ee2eb7da2a",
    "clips": "5c16ef8122f6f1798acc30eb8a30828c",
    "guess": "31b988aab6bda84867e704b9edd8a763",
}
TEXT_BAM_MD5 = "c5d3d376ca7ffc32ef5cbabcc9850804"
TEXT_PBI_MD5 = "d1d6a2f961b9fb3d29837555706c59eb"


async def check_text_input(reads: str, reports: dict[str, str], out_ext: str, reads_md5: str) -> None:
    result = await run(reads)
    await check_reports(result, reports)
    await check_bam(result, TEXT_BAM_MD5, TEXT_PBI_MD5)
    assert names(result.reads) == [f"{PAIR}.{out_ext}"], names(result.reads)
    if out_ext.endswith(".gz"):
        await assert_gunzipped_md5(result.reads[0], reads_md5, label=f"{PAIR}.{out_ext}")
    else:
        await assert_md5(result.reads[0], reads_md5, label=f"{PAIR}.{out_ext}")
    extra = [n for n in names([f async for f in result.results.walk()]) if n.endswith((".xml", ".json"))]
    assert extra == [], f"lima wrote unexpected dataset files: {extra}"


@env.task
async def test_lima_bam() -> None:
    # upstream case: lima "LIMA - Primer Removal - Input => bam"
    result = await run("bam/alz.ccs.bam")
    await check_reports(
        result,
        {
            "counts": "842c6a23ca2de504ced4538ad5111da1",
            "report": "dc073985322ae0a003ccc7e0fa4db5e6",
            "summary": "bcbcaaaca418bdeb91141c81715ca420",
            "clips": "fa03bc75bd78b2648a139fd67c69208f",
            "guess": "d3675af3ca8a908ee9e3c231668392d3",
        },
    )
    await check_bam(result, "59b04f200c309b0a60a3f182d22f6910", "851cf26eb54e4399cba5241db969dc0c")
    assert result.reads == [], f"lima wrote FASTA/FASTQ from a BAM: {names(result.reads)}"


@env.task
async def test_lima_fasta() -> None:
    # upstream case: lima "LIMA - Primer Removal - Input => fa"
    await check_text_input("fasta/alz.ccs.fasta", FASTA_REPORTS, "fasta", "d4eea419961c205ef84ef39176e30ef3")


@env.task
async def test_lima_fasta_gz() -> None:
    # upstream case: lima "LIMA - Primer Removal - Input => fa.gz"
    await check_text_input("fasta/alz.ccs.fasta.gz", FASTA_REPORTS, "fasta.gz", "d4eea419961c205ef84ef39176e30ef3")


@env.task
async def test_lima_fastq() -> None:
    # upstream case: lima "LIMA - Primer Removal - Input => fq"
    await check_text_input("fastq/alz.ccs.fastq", FASTQ_REPORTS, "fastq", "ef395f689c5566f501e300bb83d7a5f2")


@env.task
async def test_lima_fastq_gz() -> None:
    # upstream case: lima "LIMA - Primer Removal - Input => fq.gz"
    await check_text_input("fastq/alz.ccs.fastq.gz", FASTQ_REPORTS, "fastq.gz", "ef395f689c5566f501e300bb83d7a5f2")


tests = [test_lima_bam, test_lima_fasta, test_lima_fasta_gz, test_lima_fastq, test_lima_fastq_gz]
