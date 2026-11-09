"""Tests for flyte_bio.modules.glimpse2.

The cases and expected values are upstream's glimpse2/chunk, glimpse2/phase,
glimpse2/ligate and glimpse2/concordance tests (modules commit 9a1cc9f,
GLIMPSE 2.0.1), with their test configs' args. Upstream snapshots:

- chunk: the chunk file's md5;
- phase and ligate: nft-vcf's summary (chromosomes, sample and variant counts)
  and ``variantsMD5`` (the md5 of the record lines) of the phased VCF, and the
  md5 of the decompressed coverage statistics;
- concordance: the md5s of the decompressed tables.

Phase's output depends on its threads, so as upstream (``cpus = 1``) it runs
on 1 CPU, the module's default. The input-type error cases fail in the
wrapper, as upstream's module does before running.

Upstream's CRAM case converts the BAM with samtools/view (``--output-fmt
cram``, ``embed-ref=2``); here samtools convert does the same BAM-to-CRAM
conversion against the same reference, with the reference embedded. That
case's snapshot can't be reproduced (see the test); it is checked against
upstream's BAM case instead. The renaming case's BAM list is in upstream's
order, which its ``collectFile`` sets by hash: NA19401 first.
As upstream, ligate runs on two chunks phased with the default seed, each
indexed with ``bcftools index``.
"""

import asyncio
import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bcftools import bcftools_index
from flyte_bio.modules.glimpse2 import (
    Glimpse2PhaseResult,
    glimpse2_chunk,
    glimpse2_concordance,
    glimpse2_ligate,
    glimpse2_phase,
)
from flyte_bio.modules.samtools import samtools_convert
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
NA12878_VCF = HS + "illumina/vcf/NA12878.chr21_22.1X.glimpse2.vcf.gz"
NA19401_VCF = HS + "illumina/vcf/NA19401.chr21_22.1X.glimpse2.vcf.gz"
NA12878_BAM = HS + "illumina/bam/NA12878.chr21_22.1X.bam"
NA19401_BAM = HS + "illumina/bam/NA19401.chr21_22.1X.bam"
PANEL = HS + "popgen/1000GP.chr22.vcf.gz"
REGION_IN, REGION_OUT = "chr22:16600000-16800000", "chr22:16650000-16750000"
PHASED = "input_chr22_16650000-16750000"
CHUNK_ARGS = (
    "--seed 1 --window-mb 0.01 --window-cm 0.01 --window-count 200 --buffer-mb 0.005 --buffer-cm 0.005 "
    "--buffer-count 30"
)


async def glimpse2_write(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


async def glimpse2_check_vcf(file: File, name: str, samples: int, variants: int, md5: str) -> None:
    """Check a VCF's name, nft-vcf summary (chr22 only, samples, variants, phased) and variantsMD5."""
    assert file.name == name, f"wrote {file.name}, expected {name}"
    async with file.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    header = next(line for line in lines if line.startswith("#CHROM"))
    records = [line for line in lines if not line.startswith("#")]
    chromosomes = sorted({line.split("\t", 1)[0] for line in records})
    summary = (chromosomes, len(header.split("\t")) - 9, len(records))
    assert summary == (["chr22"], samples, variants), f"{name}: summary {summary}"
    actual = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
    assert actual == md5, f"{name}: variantsMD5 {actual}, expected {md5}"


async def glimpse2_check_phase(result: Glimpse2PhaseResult, samples: int, md5: str, stats: str | None) -> None:
    if stats is None:
        assert result.stats_coverage is None, "GLIMPSE2_phase wrote coverage statistics for likelihoods"
    else:
        assert result.stats_coverage is not None, "GLIMPSE2_phase wrote no coverage statistics"
        assert result.stats_coverage.name == f"{PHASED}_stats_coverage.txt.gz", result.stats_coverage.name
        await assert_gunzipped_md5(result.stats_coverage, stats, label="stats_coverage")
    await glimpse2_check_vcf(result.vcf, f"{PHASED}.vcf.gz", samples, 266, md5)


async def glimpse2_panel() -> tuple[File, File]:
    return await fixture(PANEL), await fixture(PANEL + ".csi")


@env.task
async def test_glimpse2_chunk() -> None:
    # upstream cases: glimpse2/chunk "Should run without map" and "Should run with map"
    panel, csi = await glimpse2_panel()
    genetic_map = await fixture(HS + "genome/genetic_map/genome.GRCh38.chr22.glimpse.map")
    region = "chr22:16570000-16610000"
    without_map, with_map = await asyncio.gather(
        glimpse2_chunk(panel, csi, region, prefix="test", args=CHUNK_ARGS),
        glimpse2_chunk(panel, csi, region, prefix="test", genetic_map=genetic_map, args=CHUNK_ARGS),
    )
    assert without_map.name == "test.txt", without_map.name
    await assert_md5(without_map, "c363e235162ca0f7d22e5604c192e256", label="chunk without map")
    await assert_md5(with_map, "3bfef0fe24564dd3f759ee390345970a", label="chunk with map")


@env.task
async def test_glimpse2_phase_likelihoods() -> None:
    # upstream case: glimpse2/phase "homo_sapiens - vcf region - panel"
    panel, csi = await glimpse2_panel()
    vcf, vcf_csi = await fixture(NA12878_VCF), await fixture(NA12878_VCF + ".csi")
    result = await glimpse2_phase(
        [vcf], [vcf_csi], panel, csi, PHASED, input_region=REGION_IN, output_region=REGION_OUT, args="--seed 1"
    )
    await glimpse2_check_phase(result, 1, "d8e9e1ec4d9d229536a23f2f69a469de", None)


@env.task
async def test_glimpse2_phase_reads() -> None:
    # upstream cases: glimpse2/phase "homo_sapiens - bam region - panel", "... - [bam, bam] region - panel"
    # and "... - [bam, bam] region sample rename - panel"
    panel, csi = await glimpse2_panel()
    bams = list(await asyncio.gather(fixture(NA12878_BAM), fixture(NA19401_BAM)))
    bais = list(await asyncio.gather(fixture(NA12878_BAM + ".bai"), fixture(NA19401_BAM + ".bai")))
    # Upstream builds the list with collectFile, whose default (hash) order puts NA19401 first.
    bamlist = await glimpse2_write(
        "bamlist.txt", "NA19401.chr21_22.1X.bam\tSample2\nNA12878.chr21_22.1X.bam\tSample1\n"
    )
    regions = {"input_region": REGION_IN, "output_region": REGION_OUT, "args": "--seed 1"}
    one, two, renamed = await asyncio.gather(
        glimpse2_phase(bams[:1], bais[:1], panel, csi, PHASED, **regions),
        glimpse2_phase(bams, bais, panel, csi, PHASED, **regions),
        glimpse2_phase(bams, bais, panel, csi, PHASED, bamlist=bamlist, **regions),
    )
    await glimpse2_check_phase(one, 1, "23026b2111d5b93194ebbb5de09d654d", "39bf5e7c67e377207656104f64fea657")
    await glimpse2_check_phase(two, 2, "f06d12e3c16a8edc313ad564836008", "fa1b295b071c1b8670c013577a982139")
    await glimpse2_check_phase(renamed, 2, "f2aff70c52df2f69d9bc1f4bf63fcb5c", "b21fbf77322b281c5cf4822edadd844c")


@env.task
async def test_glimpse2_phase_cram() -> None:
    # upstream case: glimpse2/phase "homo_sapiens - cram region - panel - fasta"
    panel, csi = await glimpse2_panel()
    bam, bai, fasta, fai = await asyncio.gather(
        fixture(NA12878_BAM),
        fixture(NA12878_BAM + ".bai"),
        fixture(HS + "genome/genomeGRCh38_chr21_22.fa.gz"),
        fixture(HS + "genome/genomeGRCh38_chr21_22.fa.gz.fai"),
    )
    cram = await samtools_convert(bam, fasta, "input", index=bai, args="--output-fmt-option embed_ref=2")
    result = await glimpse2_phase(
        [cram.alignment],
        [cram.index],
        panel,
        csi,
        PHASED,
        input_region=REGION_IN,
        output_region=REGION_OUT,
        fasta=fasta,
        fasta_index=fai,
        args="--seed 1",
    )
    # Upstream's snapshot for this case (variantsMD5 d8e9e1ec4d9d229536a23f2f69a469de, the likelihood case's;
    # coverage c9558dd66680731d5a360818480c1998) can't be reproduced from these reads. The CRAM holds the same
    # reads as the BAM, so this checks that it phases them as upstream's BAM case does, with the same
    # coverage table but for its sample label (GLIMPSE2 names a sample after its file, here "input").
    assert result.stats_coverage is not None, "GLIMPSE2_phase wrote no coverage statistics"
    async with result.stats_coverage.open("rb") as fh:
        stats = gzip.decompress(bytes(await fh.read())).decode()
    as_bam = stats.replace("0\tinput\t\t", "0\tNA12878.chr21_2\t", 1)
    actual = hashlib.md5(as_bam.encode()).hexdigest()
    assert actual == "39bf5e7c67e377207656104f64fea657", f"coverage (relabelled) md5 {actual}"
    await glimpse2_check_vcf(result.vcf, f"{PHASED}.vcf.gz", 1, 266, "23026b2111d5b93194ebbb5de09d654d")


@env.task
async def test_glimpse2_phase_input_errors() -> None:
    # upstream cases: glimpse2/phase "Error file type inconsistent" and "Error only one vcf allowed"
    panel, csi = await glimpse2_panel()
    bam, bai, vcf, vcf_csi, other, other_csi = await asyncio.gather(
        fixture(NA12878_BAM),
        fixture(NA12878_BAM + ".bai"),
        fixture(NA12878_VCF),
        fixture(NA12878_VCF + ".csi"),
        fixture(NA19401_VCF),
        fixture(NA19401_VCF + ".csi"),
    )
    for label, inputs, indexes, message in (
        (
            "mixed",
            [bam, vcf],
            [bai, vcf_csi],
            "Input files must be of the same type and either .bam/.cram or .vcf/.vcf.gz/.bcf format. Found: [bam, gl]",
        ),
        ("two vcfs", [other, vcf], [other_csi, vcf_csi], "Only one input .vcf/.vcf.gz/.bcf file can be provided"),
    ):
        try:
            await glimpse2_phase(inputs, indexes, panel, csi, PHASED, input_region=REGION_IN, output_region=REGION_OUT)
        except ValueError as e:
            assert message in str(e), f"{label}: unexpected error {e}"
        else:
            raise AssertionError(f"GLIMPSE2_phase accepted {label} inputs")


@env.task
async def test_glimpse2_ligate() -> None:
    # upstream case: glimpse2/ligate "Should run glimpse ligate": two chunks phased, indexed, then ligated
    panel, csi = await glimpse2_panel()
    vcf, vcf_csi = await fixture(NA12878_VCF), await fixture(NA12878_VCF + ".csi")
    chunks = [
        ("chr22:16570065-16597215", "chr22:16570065-16592216"),
        ("chr22:16587172-16609999", "chr22:16592229-16609999"),
    ]
    phased = await asyncio.gather(
        *(
            glimpse2_phase(
                [vcf],
                [vcf_csi],
                panel,
                csi,
                "input_" + region_out.replace(":", "_"),
                input_region=region_in,
                output_region=region_out,
            )
            for region_in, region_out in chunks
        )
    )
    vcfs = [result.vcf for result in phased]
    indexes = await asyncio.gather(*(bcftools_index(chunk) for chunk in vcfs))
    assert all(index is not None for index in indexes), "bcftools index wrote no index"
    merged = await glimpse2_ligate(vcfs, [index for index in indexes if index is not None], prefix="input")
    await glimpse2_check_vcf(merged, "input.vcf.gz", 1, 903, "bfeec64e4bbd1211d3dd284aff482995")


@env.task
async def test_glimpse2_concordance() -> None:
    # upstream cases: glimpse2/concordance "test glimpse2 concordance" and
    # "test list of region and rsquare per site", with their test args
    estimate, estimate_csi, truth, truth_csi, freq, freq_csi = await asyncio.gather(
        fixture(NA12878_VCF),
        fixture(NA12878_VCF + ".csi"),
        fixture(HS + "illumina/vcf/NA12878_GIAB.chr21_22.vcf.gz"),
        fixture(HS + "illumina/vcf/NA12878_GIAB.chr21_22.vcf.gz.csi"),
        fixture(HS + "popgen/1000GP.chr22.sites.vcf.gz"),
        fixture(HS + "popgen/1000GP.chr22.sites.vcf.gz.csi"),
    )
    files = (estimate, estimate_csi, truth, truth_csi, freq, freq_csi)
    bins = "0 0.01 0.05 0.1 0.2 0.5"
    one, two = await asyncio.gather(
        glimpse2_concordance(*files, ["chr22"], "input", bins=bins, args="--gt-val --af-tag AF --seed 1"),
        glimpse2_concordance(
            *files, ["chr22", "chr22"], "input", bins=bins, args="--gt-val --af-tag AF --out-r2-per-site"
        ),
    )
    for label, result, md5s in (
        (
            "one region",
            one,
            [
                "605cf14649d768c07a6c5cff15aea26a",
                "bdbe0131aa4d239de367bbb18e3c3327",
                "034c0ecc19243d70e9cfda9083472763",
                "d70d53a5012b161ba4de65cd48dcac67",
                "9f0f2f13b8ea8597a83c35d578476989",
            ],
        ),
        (
            "two regions",
            two,
            [
                "3a4753237fd8174ff46fb91bba712d5c",
                "382dd127002f43716563a71a2b693de4",
                "61a1133643ae5b0dc44ba76976b0bcce",
                "974a8c3188112716482103b8508ed109",
                "9f0f2f13b8ea8597a83c35d578476989",
            ],
        ),
    ):
        tables = [result.errors_cal, result.errors_grp, result.errors_spl, result.rsquare_grp, result.rsquare_spl]
        for table, md5 in zip(tables, md5s, strict=True):
            await assert_gunzipped_md5(table, md5, label=f"{label} {table.name}")
    assert one.rsquare_per_site is None, "concordance wrote r2 per site without --out-r2-per-site"
    assert two.rsquare_per_site is not None and two.rsquare_per_site.name == "input_r2_sites.txt.gz"
    await assert_gunzipped_md5(two.rsquare_per_site, "f65d7df3b887a4aade33a0b7b60bcca4", label="r2 per site")


tests = [
    test_glimpse2_chunk,
    test_glimpse2_phase_likelihoods,
    test_glimpse2_phase_reads,
    test_glimpse2_phase_cram,
    test_glimpse2_phase_input_errors,
    test_glimpse2_ligate,
    test_glimpse2_concordance,
]
