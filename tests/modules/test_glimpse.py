"""Tests for flyte_bio.modules.glimpse.

The cases and expected values are upstream's glimpse/chunk, glimpse/phase and
glimpse/ligate tests (modules commit 9a1cc9f, GLIMPSE 1.1.1). Upstream
snapshots the chunk file's md5, and for phased and ligated VCFs nft-vcf's
summary (chromosomes, sample and variant counts) and ``variantsMD5`` (the md5
of the record lines). As upstream, ligate runs on phase's two chunks, each
indexed with ``bcftools index``. The modules' CPUs match upstream's test CPUs
(4 for chunk and phase, 2 for ligate).
"""

import asyncio
import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bcftools import bcftools_index
from flyte_bio.modules.glimpse import glimpse_chunk, glimpse_ligate, glimpse_phase
from tests.framework import assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
SAMPLE_VCF = HS + "illumina/vcf/NA12878.chr21_22.1X.glimpse2.vcf.gz"
# Upstream's chunks: (input region, output region).
CHUNKS = [
    ("chr21:16570065-16597215", "chr21:16570065-16592216"),
    ("chr21:16587172-16609999", "chr21:16592229-16609999"),
]


async def glimpse_vcf_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode().splitlines()


async def glimpse_check_vcf(file: File, name: str, chromosome: str, samples: int, variants: int, md5: str) -> None:
    """Check a VCF's name, nft-vcf summary (chromosomes, samples, variants) and variantsMD5."""
    assert file.name == name, f"wrote {file.name}, expected {name}"
    lines = await glimpse_vcf_lines(file)
    header = next(line for line in lines if line.startswith("#CHROM"))
    records = [line for line in lines if not line.startswith("#")]
    chromosomes = sorted({line.split("\t", 1)[0] for line in records})
    summary = (chromosomes, len(header.split("\t")) - 9, len(records))
    assert summary == ([chromosome], samples, variants), f"{name}: summary {summary}"
    actual = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
    assert actual == md5, f"{name}: variantsMD5 {actual}, expected {md5}"


async def glimpse_samples_file() -> File:
    """Upstream's sampleinfos.txt: NA12878 is diploid."""
    path = Path(tempfile.mkdtemp()) / "sampleinfos.txt"
    path.write_text("NA12878 2")
    return await File.from_local(str(path))


async def glimpse_phase_chunks() -> list[File]:
    """Upstream's phase case: NA12878's chr21 likelihoods phased in its two chunks."""
    vcf, csi, panel, panel_csi, genetic_map, samples = await asyncio.gather(
        fixture(SAMPLE_VCF),
        fixture(SAMPLE_VCF + ".csi"),
        fixture(HS + "popgen/1000GP.chr21.vcf.gz"),
        fixture(HS + "popgen/1000GP.chr21.vcf.gz.csi"),
        fixture(HS + "genome/genetic_map/genome.GRCh38.chr21.glimpse.map"),
        glimpse_samples_file(),
    )
    return list(
        await asyncio.gather(
            *(
                glimpse_phase(
                    vcf,
                    csi,
                    panel,
                    panel_csi,
                    input_region,
                    output_region,
                    prefix="input_" + input_region.replace(":", "_"),
                    samples_file=samples,
                    genetic_map=genetic_map,
                )
                for input_region, output_region in CHUNKS
            )
        )
    )


@env.task
async def test_glimpse_chunk() -> None:
    # upstream case: glimpse/chunk "homo_sapiens - vcf chr 22", with its test config's window and buffer
    panel, csi = await fixture(HS + "popgen/1000GP.chr22.vcf.gz"), await fixture(HS + "popgen/1000GP.chr22.vcf.gz.csi")
    chunks = await glimpse_chunk(
        panel, csi, region="chr22", prefix="input", args="--window-size 2000000 --buffer-size 200000"
    )
    assert chunks.name == "input.txt", chunks.name
    await assert_md5(chunks, "38bc015bc7bb9f4fc0fdb08a698cfc78", label="GLIMPSE_chunk")


@env.task
async def test_glimpse_phase() -> None:
    # upstream case: glimpse/phase "test_glimpse_phase"
    first, second = await glimpse_phase_chunks()
    await glimpse_check_vcf(
        first, "input_chr21_16570065-16597215.vcf.gz", "chr21", 1, 563, "307c42e67924110b165f1fff0e5c14a"
    )
    await glimpse_check_vcf(
        second, "input_chr21_16587172-16609999.vcf.gz", "chr21", 1, 459, "b138aad7ca677f696b7635802d3be82e"
    )


@env.task
async def test_glimpse_ligate() -> None:
    # upstream case: glimpse/ligate "test_glimpse_ligate": phase's chunks, indexed, then ligated
    chunks = await glimpse_phase_chunks()
    indexes = await asyncio.gather(*(bcftools_index(chunk) for chunk in chunks))
    assert all(index is not None for index in indexes), "bcftools index wrote no index"
    merged = await glimpse_ligate(chunks, [index for index in indexes if index is not None], prefix="input")
    await glimpse_check_vcf(merged, "input.vcf.gz", "chr21", 1, 836, "bcb1deefc22f80c32d5ca5cedb56cd75")


tests = [test_glimpse_chunk, test_glimpse_phase, test_glimpse_ligate]
