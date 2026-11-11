"""Tests for flyte_bio.modules.cutesv.

The cases are upstream's cutesv tests (modules commit 9a1cc9f, cuteSV 2.1.4),
with their test config's ``-s 1`` (minimum support) for the regions case.
cuteSV's VCF header records the date and command line, so, as upstream,
each case checks the VCF's summary (chromosomes, samples, variant count,
phasing; see :func:`vcf_summary`), and the regions case checks that every
call lies in the region.
"""

import asyncio
import gzip
import tempfile
from dataclasses import dataclass
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.cutesv import cutesv
from tests.framework import env, fixture


@dataclass
class VcfSummary:
    chromosomes: list[str]  # the records' chromosomes, sorted
    sample_count: int
    variant_count: int
    phased: bool  # every record's genotypes are phased (true for no records)


async def vcf_lines(vcf: File) -> list[str]:
    """The lines of ``vcf``, plain or gzipped."""
    async with vcf.open("rb") as fh:
        data = bytes(await fh.read())
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode().splitlines()


async def vcf_summary(vcf: File) -> VcfSummary:
    """What nft-vcf's ``vcf.summary`` reports, computed from the VCF text."""
    lines = await vcf_lines(vcf)
    header = next(line for line in lines if line.startswith("#CHROM"))
    records = [line.split("\t") for line in lines if line and not line.startswith("#")]
    samples = max(len(header.split("\t")) - 9, 0)

    def phased(fields: list[str]) -> bool:
        if samples == 0:
            return False
        keys = fields[8].split(":")
        if "GT" not in keys:
            return False
        gt = keys.index("GT")
        return all("|" in s.split(":")[gt] for s in fields[9:])

    return VcfSummary(
        chromosomes=sorted({r[0] for r in records}),
        sample_count=samples,
        variant_count=len(records),
        phased=all(phased(r) for r in records),
    )


async def cutesv_text_file(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


@env.task
async def test_cutesv_sarscov2() -> None:
    # upstream case: cutesv "sarscov2 - bam - fasta"
    bam, bai, fasta = await asyncio.gather(
        fixture("genomics/sarscov2/nanopore/bam/test.sorted.bam"),
        fixture("genomics/sarscov2/nanopore/bam/test.sorted.bam.bai"),
        fixture("genomics/sarscov2/genome/genome.fasta"),
    )
    vcf = await cutesv(bam, bai, fasta, prefix="test")
    assert vcf.path.endswith("/test.vcf"), vcf.path
    summary = await vcf_summary(vcf)
    assert summary == VcfSummary([], 1, 0, True), summary


@env.task
async def test_cutesv_regions() -> None:
    # upstream case: cutesv "homo_sapiens - bam - fasta - regions"
    bam, bai, fasta, regions = await asyncio.gather(
        fixture("genomics/homo_sapiens/nanopore/bam/test.sorted.bam"),
        fixture("genomics/homo_sapiens/nanopore/bam/test.sorted.bam.bai"),
        fixture("genomics/homo_sapiens/genome/genome.fasta"),
        cutesv_text_file("regions.bed", "chr22\t18000\t21000\n"),
    )
    vcf = await cutesv(bam, bai, fasta, prefix="test", regions=regions, args="-s 1")
    summary = await vcf_summary(vcf)
    assert summary == VcfSummary(["chr22"], 1, 8, False), summary
    positions = [int(line.split("\t")[1]) for line in await vcf_lines(vcf) if not line.startswith("#")]
    assert positions and all(18000 <= p <= 21000 for p in positions), positions


tests = [test_cutesv_sarscov2, test_cutesv_regions]
