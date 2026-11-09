"""Tests for flyte_bio.modules.beagle5.

The cases and expected values are upstream's beagle5/beagle tests (modules
commit 9a1cc9f, Beagle 5.5 27Feb25.75f): the phased VCF's ``variantsMD5``
(the md5 of its record lines, as nft-vcf computes it) and the log's name.
The module's 4 CPUs match upstream's test CPUs.
"""

import asyncio
import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.beagle5 import BeagleResult, beagle5_beagle
from tests.framework import env, fixture

HS = "genomics/homo_sapiens/"
SAMPLE_VCF = HS + "illumina/vcf/NA12878.chr21_22.1X.glimpse2.vcf.gz"
PANEL = HS + "popgen/1000GP.chr22.vcf.gz"


async def beagle_variants_md5(file: File) -> str:
    async with file.open("rb") as fh:
        text = gzip.decompress(bytes(await fh.read())).decode()
    body = "".join(line for line in text.splitlines() if not line.startswith("#"))
    return format(int(hashlib.md5(body.encode()).hexdigest(), 16), "x")


async def beagle_check(result: BeagleResult, md5: str, label: str) -> None:
    names = (result.vcf.name, result.log.name)
    assert names == ("test.bglout.vcf.gz", "test.bglout.log"), f"{label}: wrote {names}"
    actual = await beagle_variants_md5(result.vcf)
    assert actual == md5, f"{label}: variantsMD5 {actual}, expected {md5}"


@env.task
async def test_beagle5_beagle() -> None:
    # upstream cases: beagle5/beagle "test-beagle5-beagle", "test-beagle5-beagle-ref" and "test-beagle5-beagle-ref-map"
    vcf, csi, panel, panel_csi, genetic_map = await asyncio.gather(
        fixture(SAMPLE_VCF),
        fixture(SAMPLE_VCF + ".csi"),
        fixture(PANEL),
        fixture(PANEL + ".csi"),
        fixture(HS + "genome/genetic_map/genome.GRCh38.chr22.plink.map"),
    )
    reference = {"reference": panel, "reference_index": panel_csi, "region": "chr22"}
    alone, with_panel, with_map = await asyncio.gather(
        beagle5_beagle(vcf, "test.bglout", vcf_index=csi),
        beagle5_beagle(vcf, "test.bglout", vcf_index=csi, **reference),
        beagle5_beagle(vcf, "test.bglout", vcf_index=csi, genetic_map=genetic_map, **reference),
    )
    await beagle_check(alone, "17649c1b1fd157f032bb238afba1df8d", "no panel")
    await beagle_check(with_panel, "9cab6e975108d5f6dbb7773b794711c2", "panel")
    await beagle_check(with_map, "9cab6e975108d5f6dbb7773b794711c2", "panel and map")


tests = [test_beagle5_beagle]
