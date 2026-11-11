"""Tests for flyte_bio.modules.svtk.

The cases are upstream's svtk/standardize and svtk/vcf2bed tests (modules
commit 9a1cc9f, svtk 0.0.20190615): the standardized VCF's md5 (decompressed,
as nf-test hashes it) and its ``##fileformat=VCFv4.2`` line, and the BED's
md5.
"""

import gzip

from flyte_bio.modules.svtk import svtk_standardize, svtk_vcf2bed
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HS = "genomics/homo_sapiens/"


@env.task
async def test_svtk_standardize_manta() -> None:
    vcf = await fixture(HS + "illumina/vcf/test.genome_21.somatic_sv.vcf.gz")
    fai = await fixture(HS + "genome/genome.fasta.fai")
    result = await svtk_standardize(vcf, "manta", prefix="test", fai=fai)
    await assert_gunzipped_md5(result, "6119075b757b3720ab5ecb484c4c8b70")
    async with result.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    assert "##fileformat=VCFv4.2" in lines, "no ##fileformat=VCFv4.2 line"


@env.task
async def test_svtk_vcf2bed() -> None:
    vcf = await fixture(HS + "illumina/vcf/NA24385_sv.vcf.gz")
    tbi = await fixture(HS + "illumina/vcf/NA24385_sv.vcf.gz.tbi")
    bed = await svtk_vcf2bed(vcf, prefix="test", tbi=tbi)
    await assert_md5(bed, "d4bff8cc165bf7213c464d90cdd326f1")


tests = [test_svtk_standardize_manta, test_svtk_vcf2bed]
