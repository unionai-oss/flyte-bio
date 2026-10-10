"""End-to-end tests for flyte_bio.pipelines.scrnaseq on the upstream test profile.

The test profile: two 10x v2 samples (Sample_Y with two lanes) against
mouse chr19 (GRCm38.p6 genome, GENCODE vM19 annotation).
"""

import asyncio
import gzip
import json

from flyte.io import File

from flyte_bio.pipelines.scrnaseq import (
    barcode_whitelist,
    prepare_index,
    prepare_reference,
    protocol,
    quantify,
    read_qc,
    samples_from_samplesheet,
)
from flyte_bio.samplesheet import read_samplesheet
from tests.framework import env, fixture, stage_samplesheet

DATA = "https://github.com/nf-core/test-datasets/raw/scrnaseq/"
SAMPLESHEET = DATA + "samplesheet-2-0.csv"
FASTA = DATA + "reference/GRCm38.p6.genome.chr19.fa"
GTF = DATA + "reference/gencode.vM19.annotation.chr19.gtf"


async def text(f: File) -> str:
    async with f.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def staged_samples():
    sheet = await stage_samplesheet(SAMPLESHEET, ("fastq_1", "fastq_2"))
    return sheet, samples_from_samplesheet(await read_samplesheet(sheet))


@env.task
async def test_prepare_reference() -> None:
    fasta, gtf = await asyncio.gather(fixture(FASTA), fixture(GTF))
    ref = await prepare_reference(fasta, gtf)
    assert ref.fasta is fasta
    assert ref.gtf.path.endswith("/GRCm38.p6.genome.chr19_genes.gtf"), ref.gtf.path
    # Every annotation line is on chr19, the genome's only sequence, so all are kept.
    original = [line for line in (await text(gtf)).splitlines() if line]
    kept = [line for line in (await text(ref.gtf)).splitlines() if line]
    assert kept == [line for line in original if line.split("\t")[0] == "chr19"], (len(kept), len(original))


@env.task
async def test_read_qc() -> None:
    _, samples = await staged_samples()
    assert [(s.id, len(s.fastq_1), s.expected_cells, s.seq_center) for s in samples] == [
        ("Sample_X", 1, 5000, "Broad Institute"),
        ("Sample_Y", 2, 5000, "CRG Barcelona"),
    ], samples
    reports = await read_qc(samples[1])
    # Two lanes of R1 + R2: four reports, numbered in lane order.
    for n in (1, 2, 3, 4):
        assert await reports.get_file(f"Sample_Y_{n}_fastqc.html") is not None, f"no Sample_Y_{n}_fastqc.html"


async def read_count(f: File) -> int:
    local = await f.download()
    with gzip.open(local, "rb") as fh:
        return sum(1 for _ in fh) // 4


@env.task
async def test_quantify() -> None:
    _, samples = await staged_samples()
    sample = samples[1]  # Sample_Y: two lanes
    fasta, gtf = await asyncio.gather(fixture(FASTA), fixture(GTF))
    proto = protocol("10XV2")
    reference, whitelist = await asyncio.gather(prepare_reference(fasta, gtf), barcode_whitelist(proto))
    assert whitelist is not None and whitelist.path.endswith("10x_V2_barcode_whitelist.txt.gz"), whitelist
    index = await prepare_index(reference)
    assert index.t2g is not None and index.t2g.path.endswith("/t2g_3col.tsv"), index.t2g
    quant = await quantify(sample, index, proto, whitelist)
    for name in ("quants_mat.mtx", "quants_mat_cols.txt", "quants_mat_rows.txt", "quants.h5ad"):
        assert await quant.quant.get_file(f"alevin/{name}") is not None, f"no alevin/{name}"
    # Both lanes were mapped: the mapper saw every read pair of the sample.
    info = json.loads(await text(await quant.map.get_file("map_info.json")))
    reads = sum(await asyncio.gather(*(read_count(f) for f in sample.fastq_1)))
    assert reads in info.values(), (reads, info)


tests = [test_prepare_reference, test_read_qc, test_quantify]
