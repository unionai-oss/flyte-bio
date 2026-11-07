"""Tests for flyte_bio.modules.pbmm2.

The cases and expected md5s are upstream's pbmm2/align tests (modules commit
9a1cc9f, pbmm2 26.2.0), with upstream's test ``args`` (``--sort --bam-index
BAI``) and its 4 CPUs. nf-test hashes the BAM and its index as they are; the
two cases' BAMs differ only in the reference name their headers record
(``genome3.fasta``, or ``genome3.fa`` for the renamed ``.fna``).
"""

import asyncio

from flyte_bio.modules.pbmm2 import pbmm2_align
from tests.framework import assert_md5, env, fixture

PACBIO_BAM = "genomics/homo_sapiens/pacbio/bam/NA03697B2_downsampled.pbmm2.repeats.bam"
GENOME = "genomics/homo_sapiens/genome/"
ARGS = "--sort --bam-index BAI"


async def check(reference: str, bam_md5: str) -> None:
    reads, fasta = await asyncio.gather(fixture(PACBIO_BAM), fixture(GENOME + reference))
    result = await pbmm2_align(reads, fasta, prefix="test", args=ARGS)
    assert result.bam.name == "test.bam", result.bam.name
    assert result.index is not None and result.index.name == "test.bam.bai", result.index
    await assert_md5(result.bam, bam_md5, label=f"{reference}: test.bam")
    await assert_md5(result.index, "d19601b196d069d281c0280bbb6dda6f", label=f"{reference}: test.bam.bai")


@env.task
async def test_pbmm2_align() -> None:
    # upstream case: pbmm2/align "pbmm2 - bam"
    await check("genome3.fasta", "21fde9abeef32bf602d7e07521d6dd58")


@env.task
async def test_pbmm2_align_fna() -> None:
    # upstream case: pbmm2/align "pbmm2 - bam - fna"
    await check("genome3.fna", "85fc597755dd2aca64e9141248fac2e0")


tests = [test_pbmm2_align, test_pbmm2_align_fna]
