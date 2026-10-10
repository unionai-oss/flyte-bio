"""End-to-end tests for flyte_bio.pipelines.scrnaseq on the upstream test profile.

The test profile: two 10x v2 samples (Sample_Y with two lanes) against
mouse chr19 (GRCm38.p6 genome, GENCODE vM19 annotation).
"""

import asyncio

from flyte.io import File

from flyte_bio.pipelines.scrnaseq import prepare_reference, read_qc, samples_from_samplesheet
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


tests = [test_prepare_reference, test_read_qc]
