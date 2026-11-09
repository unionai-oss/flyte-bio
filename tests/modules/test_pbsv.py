"""Tests for flyte_bio.modules.pbsv.

The cases and expected values are upstream's pbsv/discover and pbsv/call
tests (modules commit 9a1cc9f, pbsv 2.11.0): discover on a downsampled
PacBio BAM (nf-test hashes the ``.svsig.gz`` decompressed), and call on
upstream's signature fixture against ``genome3.fasta``, checked as upstream
does by the VCF's records (nft-vcf's ``variantsMD5``), since its header
records the run.

:func:`test_pbsv_discover_call` also chains discover into call. Upstream
has no snapshot for the chain, so it checks the VCF's structure: a pbsv
VCF over ``genome3.fasta``'s contig, holding the BAM's sample.
"""

import hashlib

from flyte.io import File

from flyte_bio.modules.pbsv import pbsv_call, pbsv_discover
from tests.framework import assert_gunzipped_md5, env, fixture

BAM = "genomics/homo_sapiens/pacbio/bam/NA03697B2_downsampled.pbmm2.repeats.bam"
SVSIG = "genomics/homo_sapiens/pacbio/svsig/NA03697B2_new.pbmm2.repeats.svsig.gz"
FASTA = "genomics/homo_sapiens/genome/genome3.fasta"


async def pbsv_lines(vcf: File) -> list[str]:
    async with vcf.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def pbsv_variants_md5(vcf: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of the record lines joined without newlines, in hex without leading zeros."""
    body = "".join(line for line in await pbsv_lines(vcf) if not line.startswith("#"))
    return format(int(hashlib.md5(body.encode()).hexdigest(), 16), "x")


@env.task
async def test_pbsv_discover_bam() -> None:
    svsig = await pbsv_discover(await fixture(BAM), prefix="test")
    assert svsig.name == "test.svsig.gz", svsig.name
    await assert_gunzipped_md5(svsig, "180f199ed523eb7e909cc1fd0a3e8cef", label=svsig.name)


@env.task
async def test_pbsv_call_svsig_fasta() -> None:
    vcf = await pbsv_call([await fixture(SVSIG)], await fixture(FASTA), prefix="test")
    assert vcf.name == "test.vcf", vcf.name
    variants = await pbsv_variants_md5(vcf)
    assert variants == "39821c95936a6f9539d3ad53e6562f99", f"test.vcf variantsMD5 {variants}"


@env.task
async def test_pbsv_discover_call() -> None:
    svsig = await pbsv_discover(await fixture(BAM), prefix="test")
    vcf = await pbsv_call([svsig], await fixture(FASTA), prefix="test_chain")
    lines = await pbsv_lines(vcf)
    assert lines[0].startswith("##fileformat=VCF"), lines[0]
    assert any("pbsv" in line for line in lines if line.startswith("##")), "no pbsv header line"
    assert any(line.startswith("##contig=<ID=chr19:45760000-45770300") for line in lines), "no genome3 contig"
    header = [line for line in lines if line.startswith("#CHROM")]
    assert header and header[0].split("\t")[9:] == ["NA03697B2"], header
    records = [line.split("\t") for line in lines if not line.startswith("#")]
    assert records, "the chain called no variants"
    assert all(r[0] == "chr19:45760000-45770300" and "SVTYPE=" in r[7] for r in records), records[:3]


tests = [
    test_pbsv_discover_bam,
    test_pbsv_call_svsig_fasta,
    test_pbsv_discover_call,
]
