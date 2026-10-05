"""End-to-end test for flyte_bio.pipelines.rnaseq on the upstream test profile.

Uses the samplesheet and references from rnaseq 3.26.0's ``test`` profile
(7 runs → 5 samples, mixed single/paired-end, two multi-run samples, a
GFP additional FASTA). Outputs aren't byte-reproducible, so this is a
run-to-green check of every per-sample and reference output.
"""

import asyncio
import csv
import io

from flyte.io import File

from flyte_bio.pipelines.rnaseq import Sample, rnaseq
from tests.framework import assert_dir_nonempty, assert_nonempty, env, fixture

DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/"
SAMPLESHEET = DATA + "samplesheet/v3.10/samplesheet_test.csv"
FASTA = DATA + "reference/genome.fasta"
GTF = DATA + "reference/genes_with_empty_tid.gtf.gz"
TRANSCRIPT_FASTA = DATA + "reference/transcriptome.fasta"
ADDITIONAL_FASTA = DATA + "reference/gfp.fa.gz"


async def load_samples(samplesheet: File) -> list[Sample]:
    """Group samplesheet rows (one per run) into Samples, keeping run order."""
    async with samplesheet.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    rows = list(csv.DictReader(io.StringIO(text)))

    async def fetch(url: str) -> File | None:
        return await fixture(url) if url else None

    r1s = await asyncio.gather(*(fetch(r["fastq_1"]) for r in rows))
    r2s = await asyncio.gather(*(fetch(r["fastq_2"]) for r in rows))

    runs: dict[str, list[tuple[File, File | None, str]]] = {}
    for row, r1, r2 in zip(rows, r1s, r2s):
        runs.setdefault(row["sample"], []).append((r1, r2, row["strandedness"]))
    return [
        Sample(
            id=sample_id,
            fastq_1=[r1 for r1, _, _ in sample_runs],
            fastq_2=[r2 for _, r2, _ in sample_runs if r2 is not None],
            strandedness=sample_runs[0][2],
        )
        for sample_id, sample_runs in runs.items()
    ]


@env.task
async def test_rnaseq_star_salmon() -> None:
    samplesheet, fasta, gtf, transcript_fasta, additional_fasta = await asyncio.gather(
        fixture(SAMPLESHEET), fixture(FASTA), fixture(GTF), fixture(TRANSCRIPT_FASTA), fixture(ADDITIONAL_FASTA)
    )
    samples = await load_samples(samplesheet)
    assert [s.id for s in samples] == [
        "WT_REP1",
        "WT_REP2",
        "RAP1_UNINDUCED_REP1",
        "RAP1_UNINDUCED_REP2",
        "RAP1_IAA_30M_REP1",
    ], [s.id for s in samples]

    result = await rnaseq(
        samples,
        fasta=fasta,
        gtf=gtf,
        transcript_fasta=transcript_fasta,
        additional_fasta=additional_fasta,
    )

    g = result.genome
    for label, f in [
        ("fasta", g.fasta),
        ("gtf", g.gtf),
        ("fai", g.fai),
        ("chrom sizes", g.chrom_sizes),
        ("gene bed", g.gene_bed),
        ("transcript fasta", g.transcript_fasta),
    ]:
        await assert_nonempty(f, label=f"genome {label}")
    await assert_dir_nonempty(g.star_index, label="star index")

    assert [r.sample for r in result.samples] == [s.id for s in samples]
    for r in result.samples:
        a = r.alignment
        for label, f in [
            ("bam", a.bam),
            ("bai", a.bai),
            ("stats", a.stats),
            ("flagstat", a.flagstat),
            ("idxstats", a.idxstats),
            ("transcriptome bam", a.transcriptome_bam),
            ("Log.final.out", a.log_final),
            ("SJ.out.tab", a.splice_junctions),
        ]:
            await assert_nonempty(f, label=f"{r.sample} {label}")
        quant = await r.salmon.get_file("quant.sf")
        assert quant is not None, f"{r.sample}: salmon results have no quant.sf"
        await assert_nonempty(quant, label=f"{r.sample} quant.sf")


tests = [test_rnaseq_star_salmon]
