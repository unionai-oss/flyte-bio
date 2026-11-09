"""Tests for flyte_bio.modules.shapeit5.

The cases and expected values are upstream's shapeit5/phasecommon and
shapeit5/ligate tests (modules commit 9a1cc9f, SHAPEIT 5.1.1). phase_common
writes BCF, so upstream snapshots its name only. ligate's snapshot is nft-vcf's
summary (chromosomes, sample and variant counts) and ``variantsMD5`` (the md5
of the record lines); as upstream, its two chunks are cut from the same VCF
with ``bcftools view`` (BCF with a ``.csi``, missing genotypes excluded).
"""

import asyncio
import gzip
import hashlib

from flyte_bio.modules.bcftools import bcftools_view
from flyte_bio.modules.shapeit5 import shapeit5_ligate, shapeit5_phasecommon
from tests.framework import env, fixture

HS = "genomics/homo_sapiens/"
PANEL = HS + "popgen/1000GP.chr22.vcf.gz"
GIAB = HS + "illumina/vcf/NA12878_GIAB.chr21_22.vcf.gz"


@env.task
async def test_shapeit5_phasecommon() -> None:
    # upstream cases: shapeit5/phasecommon "homo sapiens - vcf, no reference, no scaffold, no map" and
    # "... with map", with the test config's --seed 1
    panel, csi, genetic_map = await asyncio.gather(
        fixture(PANEL), fixture(PANEL + ".csi"), fixture(HS + "genome/genetic_map/genome.GRCh38.chr22.glimpse.map")
    )
    without_map, with_map = await asyncio.gather(
        shapeit5_phasecommon(panel, csi, "chr22", "input", args="--seed 1"),
        shapeit5_phasecommon(panel, csi, "chr22", "input", genetic_map=genetic_map, args="--seed 1"),
    )
    assert without_map.name == "input.bcf" and with_map.name == "input.bcf", (without_map.name, with_map.name)


@env.task
async def test_shapeit5_ligate() -> None:
    # upstream case: shapeit5/ligate "homo sapiens - map, [vcf], [tbi]"
    vcf, csi = await fixture(GIAB), await fixture(GIAB + ".csi")
    regions = ["chr22:16570000-16600000", "chr22:16580000-16610000"]
    chunks = await asyncio.gather(
        *(
            bcftools_view(
                vcf,
                f"NA12878_{region}",
                index=csi,
                exclude='GT="./."||GT="."',
                args=f"--regions {region} -Ob --write-index=csi",
            )
            for region in regions
        )
    )
    assert all(chunk.index is not None for chunk in chunks), "bcftools view wrote no index"
    merged = await shapeit5_ligate(
        [chunk.vcf for chunk in chunks], [chunk.index for chunk in chunks if chunk.index is not None], "NA12878"
    )
    assert merged.name == "NA12878.vcf.gz", merged.name
    async with merged.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    header = next(line for line in lines if line.startswith("#CHROM"))
    records = [line for line in lines if not line.startswith("#")]
    summary = (sorted({line.split("\t", 1)[0] for line in records}), len(header.split("\t")) - 9, len(records))
    assert summary == (["chr22"], 1, 81), f"summary {summary}"
    variants = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
    assert variants == "e0d8b01df09f70f7ccd010f42a7e07c5", f"variantsMD5 {variants}"


tests = [test_shapeit5_phasecommon, test_shapeit5_ligate]
