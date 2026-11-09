"""Tests for flyte_bio.modules.stitch.

The cases and expected values are upstream's stitch tests (modules commit
9a1cc9f, STITCH 1.7.3), on 4 CPUs as upstream's tests ran. Upstream snapshots
the names of the input, RData and plot files, and the VCF's sample names and
nft-vcf summary (chromosomes, sample and variant counts); with a seed it also
snapshots the VCF's ``variantsMD5`` (the md5 of the record lines), or the
BGEN's name. The two-stage case generates the input files, then imputes from
them, as upstream.
"""

import asyncio
import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.stitch import StitchResult, stitch
from tests.framework import env, fixture

HS = "genomics/homo_sapiens/"
CHR21 = HS + "genome/chr21/sequence/"
READ_NAMES = ["test.paired_end.recalibrated.sorted", "test2.paired_end.recalibrated.sorted"]
REGION = "chr21.8522360.8522670"
INPUT_FILES = [f"sample.{i}.input.{REGION}.RData" for i in (1, 2)]
RDATA_FILES = [
    f"EM.all.{REGION}.RData",
    f"EM.all.{REGION}.withBuffer.RData",
    f"end.{REGION}.RData",
    f"endEM.{REGION}.RData",
    f"sampleNames.{REGION}.RData",
    f"start.{REGION}.RData",
    f"startEM.{REGION}.RData",
]
PLOT_FILES = [
    f"alphaMat.{REGION}.all.s.1.png",
    f"alphaMat.{REGION}.normalized.s.1.png",
    f"hapSum.{REGION}.s.1.png",
    f"hapSum_log.{REGION}.s.1.png",
    f"metricsForPostImputationQC.{REGION}.sample.jpg",
    f"metricsForPostImputationQCChromosomeWide.{REGION}.sample.jpg",
    f"r2.{REGION}.goodonly.jpg",
]
SEEDED_VARIANTS = "e41c51d08be2cdacd81fc8800cffb156"


async def stitch_write(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


async def stitch_params() -> dict:
    """Upstream's region parameters: 10 dbSNP sites of chr21, K=2, nGen=1, no buffer."""
    posfile, genetic_map = await asyncio.gather(
        fixture(CHR21 + "dbsnp_138.hg38.first_10_biallelic_sites.tsv"),
        fixture(CHR21 + "genome.GRCh38.stitch.chr21.txt.gz"),
    )
    return {
        "posfile": posfile,
        "genetic_map": genetic_map,
        "chromosome": "chr21",
        "start": 8522360,
        "end": 8522670,
        "buffer": 0,
        "k": 2,
        "ngen": 1,
        "prefix": "test_reads",
    }


async def stitch_reads(extension: str) -> dict:
    """The two samples' reads (``cram`` or ``bam``), their list and the chr21 reference."""
    index = "crai" if extension == "cram" else "bai"
    paths = [f"{HS}illumina/{extension}/{name}.{extension}" for name in READ_NAMES]
    reads = list(await asyncio.gather(*(fixture(p) for p in paths)))
    indexes = list(await asyncio.gather(*(fixture(f"{p}.{index}") for p in paths)))
    reads_list = await stitch_write("cramlist.txt", "".join(f"{name}.{extension}\n" for name in sorted(READ_NAMES)))
    fasta, fai = await fixture(CHR21 + "genome.fasta"), await fixture(CHR21 + "genome.fasta.fai")
    return {"reads": reads, "indexes": indexes, "reads_list": reads_list, "fasta": fasta, "fasta_index": fai}


async def stitch_check(
    result: StitchResult, samples: list[str], variants: str | None, rdata: list[str] = RDATA_FILES
) -> None:
    """Check the input, RData and plot files' names, and the VCF's samples, summary and (with a seed) variantsMD5."""
    for label, files, names in (("input", result.input, INPUT_FILES), ("RData", result.rdata, rdata)):
        actual = [f.name for f in files]
        assert actual == names, f"{label}: {actual}"
    actual_plots = [f.name for f in result.plots]
    assert actual_plots == PLOT_FILES, f"plots: {actual_plots}"
    assert result.vcf is not None and result.vcf.name == "test_reads.vcf.gz", f"vcf: {result.vcf}"
    async with result.vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    header = next(line for line in lines if line.startswith("#CHROM"))
    records = [line for line in lines if not line.startswith("#")]
    assert sorted(header.split("\t")[9:]) == samples, f"samples: {header.split(chr(9))[9:]}"
    summary = (sorted({line.split("\t", 1)[0] for line in records}), len(records))
    assert summary == (["chr21"], 10), f"summary {summary}"
    if variants is not None:
        actual = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
        assert actual == variants, f"variantsMD5 {actual}, expected {variants}"


@env.task
async def test_stitch_no_seed() -> None:
    # upstream case: stitch "test_no_seed" (unseeded, so upstream snapshots no variantsMD5)
    result = await stitch(**await stitch_reads("cram"), **await stitch_params())
    await stitch_check(result, ["normal", "tumour"], None)


@env.task
async def test_stitch_with_seed() -> None:
    # upstream case: stitch "test_with_seed": BGEN output
    result = await stitch(**await stitch_reads("cram"), **await stitch_params(), seed=1, args="--output_format bgen")
    for label, files, names in (
        ("input", result.input, INPUT_FILES),
        ("RData", result.rdata, RDATA_FILES),
        ("plots", result.plots, PLOT_FILES),
    ):
        actual = [f.name for f in files]
        assert actual == names, f"{label}: {actual}"
    assert result.vcf is None and result.bgen is not None and result.bgen.name == "test_reads.bgen", result.bgen


@env.task
async def test_stitch_two_stage_imputation() -> None:
    # upstream case: stitch "test_two_stage_imputation": generate the input files, then impute from them
    params = await stitch_params()
    first = await stitch(**await stitch_reads("cram"), **params, seed=1, args="--generateInputOnly TRUE")
    assert first.vcf is None and [f.name for f in first.input] == INPUT_FILES, [f.name for f in first.input]
    second = await stitch(
        **params,
        input_dir=first.input_dir,
        rdata_dir=first.rdata_dir,
        seed=1,
        args="--regenerateInputWithDefaultValues TRUE",
    )
    # The second stage's RData starts as a copy of the first's, which includes the positions STITCH saves
    # when only generating input (pos.<region>.RData). Upstream's snapshot lacks it: its test sets the
    # second stage's args as global params, which (its config reading them lazily) its first stage gets too,
    # so that stage runs in full and saves no positions file.
    rdata = sorted([*RDATA_FILES, f"pos.{REGION}.RData"])
    await stitch_check(second, ["normal", "tumour"], SEEDED_VARIANTS, rdata=rdata)


@env.task
async def test_stitch_with_bam() -> None:
    # upstream cases: stitch "test_with_bam" and "test_with_bam and renaming"
    reads, params = await stitch_reads("bam"), await stitch_params()
    names = await stitch_write("samplename.txt", "MySample1\nMySample2\n")
    plain, renamed = await asyncio.gather(
        stitch(**reads, **params, seed=1), stitch(**reads, **params, sample_names=names, seed=1)
    )
    await stitch_check(plain, ["normal", "tumour"], SEEDED_VARIANTS)
    await stitch_check(renamed, ["MySample1", "MySample2"], SEEDED_VARIANTS)


tests = [test_stitch_no_seed, test_stitch_with_seed, test_stitch_two_stage_imputation, test_stitch_with_bam]
