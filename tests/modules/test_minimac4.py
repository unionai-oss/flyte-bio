"""Tests for flyte_bio.modules.minimac4.

The cases and expected values are upstream's minimac4/compressref and
minimac4/impute tests (modules commit 9a1cc9f, Minimac4 4.1.6). Upstream
snapshots the compressed panel's name, and the imputed VCF's and sites file's
``variantsMD5`` (the md5 of the record lines, as nft-vcf computes it). As
upstream, impute runs on the panel compressref writes.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.minimac4 import minimac4_compressref, minimac4_impute
from tests.framework import env, fixture

HS = "genomics/homo_sapiens/"
PANEL = HS + "popgen/1000GP.chr22.vcf.gz"
TARGET = HS + "illumina/vcf/NA12878.chr21_22.1X.glimpse2.vcf.gz"
SITES_MD5 = "286e107e0d5f91f5c4719030186574bc"
IMPUTED_MD5 = "a0615a8a5ce53f496ac691c5fc82a44"


async def minimac4_variants_md5(file: File) -> str:
    async with file.open("rb") as fh:
        raw = bytes(await fh.read())
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    body = "".join(line for line in raw.decode().splitlines() if not line.startswith("#"))
    return format(int(hashlib.md5(body.encode()).hexdigest(), 16), "x")


async def minimac4_compressed_panel() -> File:
    return await minimac4_compressref(await fixture(PANEL), "input", reference_index=await fixture(PANEL + ".csi"))


@env.task
async def test_minimac4_compressref() -> None:
    # upstream case: minimac4/compressref "homo sapiens - vcf"
    msav = await minimac4_compressed_panel()
    assert msav.name == "input.msav", msav.name


@env.task
async def test_minimac4_impute() -> None:
    # upstream cases: minimac4/impute "homo sapiens - vcf, ref, sites, map" (--chunk 1000000)
    # and "homo sapiens - vcf, ref, [], []" (--output-format vcf)
    msav, target, csi, genetic_map = await asyncio.gather(
        minimac4_compressed_panel(),
        fixture(TARGET),
        fixture(TARGET + ".csi"),
        fixture(HS + "genome/genetic_map/genome.GRCh38.chr22.minimac.map"),
    )
    with_map, plain = await asyncio.gather(
        minimac4_impute(
            target,
            msav,
            "NA12878",
            target_index=csi,
            genetic_map=genetic_map,
            region="chr22",
            write_sites=True,
            args="--chunk 1000000",
        ),
        minimac4_impute(target, msav, "NA12878", target_index=csi, write_sites=True, args="--output-format vcf"),
    )
    for label, result, name in (("map", with_map, "NA12878.vcf.gz"), ("plain", plain, "NA12878.vcf")):
        assert result.vcf.name == name, f"{label}: wrote {result.vcf.name}"
        assert result.sites is not None and result.sites.name == "NA12878.sites.vcf.gz", f"{label}: no sites"
        imputed, sites = await minimac4_variants_md5(result.vcf), await minimac4_variants_md5(result.sites)
        assert imputed == IMPUTED_MD5, f"{label}: imputed variantsMD5 {imputed}"
        assert sites == SITES_MD5, f"{label}: sites variantsMD5 {sites}"


tests = [test_minimac4_compressref, test_minimac4_impute]
