"""Tests for flyte_bio.modules.hificnv.

The cases and expected md5s are upstream's hificnv tests (modules commit
9a1cc9f, HiFiCNV 1.0.1), with their config's ``--cov-regex "."`` (unquoted
here, as ``args`` isn't parsed by a shell): a downsampled PacBio BAM against
``genome3.fasta`` with no optional inputs, with the sample's VCF as ``maf``,
with an exclude BED too, and with an expected copy-number BED too. The BEDs
are upstream's one-line files (no trailing newline). nf-test hashes the
``.vcf.gz`` decompressed. As upstream, all four cases give the same tracks
and calls.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.hificnv import HificnvResult, hificnv
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

PACBIO = "genomics/homo_sapiens/pacbio/"
BAM = PACBIO + "bam/NA03697B2_downsampled.pbmm2.repeats.bam"
MAF = PACBIO + "vcf/NA03697B2_new.pbmm2.repeats.vcf.gz"
FASTA = "genomics/homo_sapiens/genome/genome3.fasta"
ARGS = "--cov-regex ."
EXCLUDE_BED = "chr19:45760000-45770300\t45760000\t45770300\tchr19_bogus\t2"
EXPECTED_CN_BED = "chr19:45760000-45770300\t45761000\t45762000\tchr19_bogus_exp\t0"


async def hificnv_bed(name: str, content: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(content)
    return await File.from_local(str(path))


async def run(prefix: str, with_maf: bool, exclude: bool = False, expected_cn: bool = False) -> HificnvResult:
    return await hificnv(
        await fixture(BAM),
        await fixture(BAM + ".bai"),
        await fixture(FASTA),
        maf=await fixture(MAF) if with_maf else None,
        exclude=await hificnv_bed("exclude.bed", EXCLUDE_BED) if exclude else None,
        expected_cn=await hificnv_bed("expected_cn.bed", EXPECTED_CN_BED) if expected_cn else None,
        prefix=prefix,
        args=ARGS,
    )


async def check(result: HificnvResult, prefix: str, with_maf: bool) -> None:
    stem = f"{prefix}.NA03697B2"
    names = [result.copynum.name, result.depth.name, result.vcf.name]
    assert names == [f"{stem}.copynum.bedgraph", f"{stem}.depth.bw", f"{stem}.vcf.gz"], names
    await assert_md5(result.copynum, "93cf41145d68bc4bf0bf679f0518505b", label=result.copynum.name)
    await assert_md5(result.depth, "5bd9a8bc063efb2020ded50ba38b4348", label=result.depth.name)
    await assert_gunzipped_md5(result.vcf, "0a653e775977d9742251b3cac5161958", label=result.vcf.name)
    if with_maf:
        assert result.maf is not None and result.maf.name == f"{stem}.maf.bw", result.maf
        await assert_md5(result.maf, "594ad937321d34420fde12ce06719efd", label=result.maf.name)
    else:
        assert result.maf is None, result.maf


@env.task
async def test_hificnv_no_optional_inputs() -> None:
    await check(await run("test", with_maf=False), "test", with_maf=False)


@env.task
async def test_hificnv_with_maf_input() -> None:
    await check(await run("test_with_maf", with_maf=True), "test_with_maf", with_maf=True)


@env.task
async def test_hificnv_with_exclude_bed() -> None:
    result = await run("test_with_exclude", with_maf=True, exclude=True)
    await check(result, "test_with_exclude", with_maf=True)


@env.task
async def test_hificnv_with_expected_cn_bed() -> None:
    result = await run("test_with_expected_cn", with_maf=True, exclude=True, expected_cn=True)
    await check(result, "test_with_expected_cn", with_maf=True)


tests = [
    test_hificnv_no_optional_inputs,
    test_hificnv_with_maf_input,
    test_hificnv_with_exclude_bed,
    test_hificnv_with_expected_cn_bed,
]
