"""Tests for flyte_bio.modules.rtgtools.

The cases are upstream's rtgtools/format tests (modules commit 9a1cc9f, RTG
Tools 3.13), whose snapshots list the files of the SDF directory: sarscov2's
genome FASTA, its paired FASTQs (``--format fastq``) and its single-end BAM
with a read group (``--format sam-se``).

As upstream, the FASTQ case passes both files but leaves ``single_end`` at its
default, so only ``test_1.fastq.gz`` is formatted.
"""

import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.rtgtools import rtgtools_format
from tests.framework import env, fixture

COMMON = [
    "done",
    "format.log",
    "mainIndex",
    "nameIndex0",
    "namedata0",
    "namepointer0",
    "progress",
    "seqdata0",
    "seqpointer0",
    "sequenceIndex0",
    "summary.txt",
]
SUFFIX = ["suffixIndex0", "suffixdata0", "suffixpointer0"]


async def sdf_files(sdf: Dir) -> list[str]:
    return sorted(f.name for f in await sdf.list_files())


async def check_sdf(sdf: Dir, expected: list[str]) -> None:
    assert sdf.path.rstrip("/").endswith("/test.sdf"), sdf.path
    names = await sdf_files(sdf)
    assert names == sorted(expected), f"SDF files {names}"


@env.task
async def test_rtgtools_format_fasta() -> None:
    sdf = await rtgtools_format(await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test")
    await check_sdf(sdf, COMMON + SUFFIX)


@env.task
async def test_rtgtools_format_fastqs() -> None:
    sdf = await rtgtools_format(
        await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz"),
        await fixture("genomics/sarscov2/illumina/fastq/test_2.fastq.gz"),
        prefix="test",
        args="--format fastq",
    )
    await check_sdf(sdf, [*COMMON, "qualitydata0", *SUFFIX])


@env.task
async def test_rtgtools_format_bam_rg() -> None:
    rg = Path(tempfile.mkdtemp()) / "sam_rg.txt"
    rg.write_text("@RG\tID:READGROUP1\tSM:SAMPLE\tPL:ILLUMINA")
    sdf = await rtgtools_format(
        await fixture("genomics/sarscov2/illumina/bam/test.single_end.bam"),
        sam_rg=await File.from_local(str(rg)),
        prefix="test",
        args="--format sam-se",
    )
    await check_sdf(sdf, [*COMMON, "qualitydata0"])


tests = [test_rtgtools_format_fasta, test_rtgtools_format_fastqs, test_rtgtools_format_bam_rg]
