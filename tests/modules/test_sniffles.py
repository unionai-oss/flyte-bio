"""Tests for flyte_bio.modules.sniffles.

The cases and expected lines are upstream's sniffles tests (modules commit
9a1cc9f, Sniffles 2.8.1), with their test config: ``--minsupport 1 --no-qc
--mapq 0`` and prefix ``<id>`` for calling, prefix ``<id>_genotyped`` and no
options for genotyping. Sniffles's VCF records its command line and date, so
as upstream the calling cases compare the VCF's lines 1-2 and 5-63 and the
genotyping cases every line but ``##genotype…``; the expected lines are kept
as the md5 of their newline-joined text. The ``.snf`` and ``.tbi`` are
checked by name, as upstream.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.sniffles import SnifflesResult, sniffles
from tests.framework import env, fixture

HUMAN = "genomics/homo_sapiens/"
GENOTYPE = HUMAN + "nanopore/bam/HG002_ont_sv_genotype/"
CALL_ARGS = "--minsupport 1 --no-qc --mapq 0"
CALL_HEAD = ["##fileformat=VCFv4.2", "##source=Sniffles2_2.8.1"]
CALL_BODY_MD5 = "ffdbcbfa7822e2472e463e261670696f"  # lines 5-63
GENOTYPED_MD5 = "c281f327e0b88fc09123713ad1e0ca75"  # all 244 lines but ##genotype…


async def sniffles_vcf_lines(vcf: File) -> list[str]:
    async with vcf.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode().splitlines()


def lines_md5(lines: list[str]) -> str:
    return hashlib.md5("\n".join(lines).encode()).hexdigest()


async def sniffles_call(vcf_output: bool, snf_output: bool, tandem: bool = False) -> SnifflesResult:
    bam, bai, fasta, bed = await asyncio.gather(
        fixture(HUMAN + "nanopore/bam/test.sorted.bam"),
        fixture(HUMAN + "nanopore/bam/test.sorted.bam.bai"),
        fixture(HUMAN + "genome/genome.fasta"),
        fixture(HUMAN + "genome/genome.bed"),
    )
    return await sniffles(
        bam,
        prefix="test",
        index=bai,
        fasta=fasta,
        tandem_file=bed if tandem else None,
        vcf_output=vcf_output,
        snf_output=snf_output,
        args=CALL_ARGS,
    )


async def check_calls(result: SnifflesResult) -> None:
    assert result.vcf is not None and result.tbi is not None, result
    assert result.tbi.path.endswith("/test.vcf.gz.tbi"), result.tbi.path
    lines = await sniffles_vcf_lines(result.vcf)
    assert lines[0:2] == CALL_HEAD, lines[0:2]
    assert len(lines[4:63]) == 59 and lines_md5(lines[4:63]) == CALL_BODY_MD5, lines[4:63]


def check_snf(result: SnifflesResult) -> None:
    assert result.snf is not None and result.snf.path.endswith("/test.snf"), result.snf


@env.task
async def test_sniffles_vcf() -> None:
    # upstream case: sniffles "homo_sapiens - [bam, bai], fasta, [], true, false"
    result = await sniffles_call(vcf_output=True, snf_output=False)
    assert result.snf is None, result.snf
    await check_calls(result)


@env.task
async def test_sniffles_snf() -> None:
    # upstream case: sniffles "homo_sapiens - [bam, bai], fasta, [], false, true"
    result = await sniffles_call(vcf_output=False, snf_output=True)
    assert result.vcf is None and result.tbi is None, result
    check_snf(result)


@env.task
async def test_sniffles_vcf_snf() -> None:
    # upstream case: sniffles "homo_sapiens - [bam, bai], fasta, [], true, true"
    result = await sniffles_call(vcf_output=True, snf_output=True)
    check_snf(result)
    await check_calls(result)


@env.task
async def test_sniffles_tandem_repeats() -> None:
    # upstream case: sniffles "homo_sapiens - [bam, bai], fasta, bed, true, true"
    result = await sniffles_call(vcf_output=True, snf_output=True, tandem=True)
    check_snf(result)
    await check_calls(result)


async def check_genotyped(sites: File) -> None:
    bam, bai = await asyncio.gather(
        fixture(GENOTYPE + "HG002_ont_sv_genotype.bam"), fixture(GENOTYPE + "HG002_ont_sv_genotype.bam.bai")
    )
    result = await sniffles(bam, prefix="hg002_genotyped", index=bai, genotype_vcf=sites)
    assert result.vcf is not None and result.tbi is not None and result.snf is None, result
    assert result.tbi.path.endswith("/hg002_genotyped.vcf.gz.tbi"), result.tbi.path
    lines = [line for line in await sniffles_vcf_lines(result.vcf) if not line.startswith("##genotype")]
    assert len(lines) == 244 and lines_md5(lines) == GENOTYPED_MD5, lines


@env.task
async def test_sniffles_genotype_vcf_gz() -> None:
    # upstream case: sniffles "homo_sapiens - nanopore [bam, bai, vcf.gz], [], [], true, false - genotype"
    await check_genotyped(await fixture(GENOTYPE + "HG002_ont_sv_sites.vcf.gz"))


@env.task
async def test_sniffles_genotype_vcf() -> None:
    # upstream case: sniffles "homo_sapiens - nanopore [bam, bai, vcf], [], [], true, false - genotype"
    await check_genotyped(await gunzip(await fixture(GENOTYPE + "HG002_ont_sv_sites.vcf.gz")))


tests = [
    test_sniffles_vcf,
    test_sniffles_snf,
    test_sniffles_vcf_snf,
    test_sniffles_tandem_repeats,
    test_sniffles_genotype_vcf_gz,
    test_sniffles_genotype_vcf,
]
