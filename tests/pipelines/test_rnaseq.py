"""End-to-end tests for flyte_bio.pipelines.rnaseq on the upstream test profile.

Uses the samplesheet and references from rnaseq 3.26.0's ``test`` profile
(7 runs → 5 samples, mixed single/paired-end, two multi-run samples, a
GFP additional FASTA). Outputs aren't byte-reproducible, so these are
run-to-green checks of every per-sample, merged and reference output: the
default path, and upstream's UMI test (``--umi_dedup_tool 'umitools'``).
"""

import asyncio
import csv
import io
import tempfile
from dataclasses import fields
from pathlib import Path

import flyte
from flyte.io import File

from flyte_bio.modules.rseqc import DEFAULT_MODULES as RSEQC_MODULES
from flyte_bio.pipelines.rnaseq import Sample, rnaseq, samples_from_samplesheet
from flyte_bio.samplesheet import read_samplesheet
from tests.framework import assert_dir_nonempty, assert_nonempty, env, fixture, fixture_dir

DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/"
SAMPLESHEET = DATA + "samplesheet/v3.10/samplesheet_test.csv"
FASTA = DATA + "reference/genome.fasta"
GTF = DATA + "reference/genes_with_empty_tid.gtf.gz"
TRANSCRIPT_FASTA = DATA + "reference/transcriptome.fasta"
ADDITIONAL_FASTA = DATA + "reference/gfp.fa.gz"
SALMON_INDEX = DATA + "reference/salmon.tar.gz"
BBSPLIT_FASTA_LIST = DATA + "reference/bbsplit_fasta_list.txt"


async def load_bbsplit_refs(fasta_list: File) -> dict[str, File]:
    """Read upstream's ``name,fasta_url`` BBSplit list and fetch each FASTA."""
    async with fasta_list.open("rb") as fh:
        rows = [r for r in csv.reader(io.StringIO(bytes(await fh.read()).decode())) if r]
    files = await asyncio.gather(*(fixture(url) for _, url in rows))
    return {name: f for (name, _), f in zip(rows, files)}


async def stage_samplesheet(url: str) -> tuple[File, list[Sample]]:
    """Upstream's samplesheet with every FASTQ swapped for its cached fixture copy.

    Returns the staged sheet (what rnaseq takes) and the Samples it reads from it.
    """
    rows = await read_samplesheet(await fixture(url))

    async def fetch(fastq: str) -> str:
        return (await fixture(fastq)).path if fastq else ""

    r1s = await asyncio.gather(*(fetch(r["fastq_1"]) for r in rows))
    r2s = await asyncio.gather(*(fetch(r["fastq_2"]) for r in rows))
    out = Path(tempfile.mkdtemp(prefix="samplesheet_")) / "samplesheet.csv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        for row, r1, r2 in zip(rows, r1s, r2s):
            writer.writerow([row["sample"], r1, r2, row["strandedness"]])
    staged = await File.from_local(str(out))
    return staged, samples_from_samplesheet(await read_samplesheet(staged))


@env.task
async def test_rnaseq_star_salmon() -> None:
    (samplesheet, samples), fasta, gtf, transcript_fasta, additional_fasta = await asyncio.gather(
        stage_samplesheet(SAMPLESHEET),
        fixture(FASTA),
        fixture(GTF),
        fixture(TRANSCRIPT_FASTA),
        fixture(ADDITIONAL_FASTA),
    )
    salmon_index = await fixture_dir(SALMON_INDEX, "salmon")
    bbsplit_refs = await load_bbsplit_refs(await fixture(BBSPLIT_FASTA_LIST))
    assert sorted(bbsplit_refs) == ["human", "sarscov2"], sorted(bbsplit_refs)
    assert [s.id for s in samples] == [
        "WT_REP1",
        "WT_REP2",
        "RAP1_UNINDUCED_REP1",
        "RAP1_UNINDUCED_REP2",
        "RAP1_IAA_30M_REP1",
    ], [s.id for s in samples]

    result = await rnaseq(
        samplesheet,
        fasta=fasta,
        gtf=gtf,
        transcript_fasta=transcript_fasta,
        additional_fasta=additional_fasta,
        salmon_index_dir=salmon_index,
        bbsplit_fasta_list=bbsplit_refs,
        skip_bbsplit=False,  # the upstream test profile turns BBSplit on
        publish_results=True,
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

    assert result.failed_trimming == {}, result.failed_trimming
    assert [r.sample for r in result.samples] == [s.id for s in samples]
    # Strandedness: only WT_REP1 is 'auto'; the GSE110004 libraries are
    # reverse-stranded (every other sample is declared 'reverse').
    for r in result.samples:
        if r.sample == "WT_REP1":
            assert r.strandedness_analysis is not None
            await assert_nonempty(r.strandedness_analysis.lib_format_counts, label="WT_REP1 lib_format_counts")
            assert r.strandedness_analysis.inferred == "reverse", r.strandedness_analysis
        else:
            assert r.strandedness_analysis is None
        assert r.strandedness == "reverse", (r.sample, r.strandedness)

    for r in result.samples:
        assert r.mapping_passed and r.percent_mapped >= 5, (r.sample, r.percent_mapped)

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
        pre = r.preprocessing
        assert pre.raw_fastqc is not None and pre.trimming is not None
        mates = ["_1", "_2"] if pre.reads_2 is not None else [""]
        for mate in mates:
            for kind in ("raw", "trimmed"):
                qc = pre.raw_fastqc if kind == "raw" else pre.trimming.results
                stem = f"{r.sample}_raw{mate}" if kind == "raw" else (
                    f"{r.sample}_trimmed{mate}_val{mate}" if mate else f"{r.sample}_trimmed_trimmed"
                )
                assert await qc.get_file(f"{stem}_fastqc.html") is not None, f"{r.sample}: no {stem}_fastqc.html"
        for report in pre.trimming.reports:
            await assert_nonempty(report, label=f"{r.sample} trimming report")
        assert pre.reads_after_trimming is not None and pre.reads_after_trimming >= 10000, pre.reads_after_trimming
        assert sorted(pre.lint) == ["bbsplit", "raw", "trimmed"], sorted(pre.lint)
        for stage, log in pre.lint.items():
            await assert_nonempty(log, label=f"{r.sample} fq lint {stage}")
        assert pre.bbsplit is not None
        await assert_nonempty(pre.bbsplit.stats, label=f"{r.sample} bbsplit stats")
        await assert_nonempty(pre.bbsplit.reads_1, label=f"{r.sample} bbsplit primary R1")
        assert (pre.bbsplit.reads_2 is None) == (pre.reads_2 is None)
        assert pre.reads_1 is pre.bbsplit.reads_1  # alignment uses the primary-genome reads
        md = r.markduplicates
        assert md is not None
        for label, f in [("bam", md.bam), ("bai", md.bai), ("metrics", md.metrics), ("flagstat", md.flagstat)]:
            await assert_nonempty(f, label=f"{r.sample} markdup {label}")
        assert r.dupradar is not None
        for label, f in [("dupMatrix", r.dupradar.dup_matrix), ("intercept mqc", r.dupradar.intercept_mqc)]:
            await assert_nonempty(f, label=f"{r.sample} dupradar {label}")
        async with r.dupradar.intercept_mqc.open("rb") as fh:
            last = bytes(await fh.read()).decode().strip().splitlines()[-1]
        assert last.split()[0] == r.sample, f"dupradar mqc sample name {last!r}"
        assert r.qualimap is not None
        for name in ("qualimapReport.html", "rnaseq_qc_results.txt"):
            assert await r.qualimap.get_file(name) is not None, f"{r.sample}: qualimap wrote no {name}"
        assert sorted(r.rseqc) == sorted(RSEQC_MODULES), sorted(r.rseqc)
        infer = await r.rseqc["infer_experiment"].get_file(f"{r.sample}.infer_experiment.txt")
        assert infer is not None, f"{r.sample}: no infer_experiment.txt"
        await assert_nonempty(infer, label=f"{r.sample} infer_experiment")
        assert r.biotype_counts is not None and r.biotype_qc is not None
        await assert_nonempty(r.biotype_counts.counts, label=f"{r.sample} biotype featureCounts")
        await assert_nonempty(r.biotype_qc.counts_mqc, label=f"{r.sample} biotype_counts_mqc")
        assert r.stringtie is not None and r.stringtie.coverage_gtf is not None
        for label, f in [("transcripts.gtf", r.stringtie.transcripts_gtf), ("abundance", r.stringtie.abundance)]:
            await assert_nonempty(f, label=f"{r.sample} stringtie {label}")
        # Every sample here is stranded (reverse), so it gets per-strand + combined tracks.
        assert sorted(r.bigwig) == ["combined", "forward", "reverse"], sorted(r.bigwig)
        for strand, f in r.bigwig.items():
            await assert_nonempty(f, label=f"{r.sample} {strand} bigWig")
        quant = await r.salmon.get_file("quant.sf")
        assert quant is not None, f"{r.sample}: salmon results have no quant.sf"
        await assert_nonempty(quant, label=f"{r.sample} quant.sf")

    m = result.salmon
    await assert_nonempty(m.tx2gene, label="salmon.merged tx2gene")
    for field in fields(m.tximport):
        await assert_nonempty(getattr(m.tximport, field.name), label=f"salmon.merged {field.name}")
    await assert_nonempty(m.gene_rds, label="gene SummarizedExperiment")
    await assert_nonempty(m.transcript_rds, label="transcript SummarizedExperiment")

    assert result.deseq2_qc is not None
    await assert_nonempty(result.deseq2_qc.pca_multiqc, label="deseq2 pca mqc")
    await assert_nonempty(result.deseq2_qc.dists_multiqc, label="deseq2 dists mqc")

    # MultiQC: one report covering every sample, with upstream's strandedness checks.
    assert result.multiqc is not None
    await assert_nonempty(result.multiqc.report, label="multiqc_report.html")
    data_files = {f.path.rsplit("/", 1)[-1]: f async for f in result.multiqc.results.walk() if "_data/" in f.path}
    assert "multiqc_general_stats.txt" in data_files, sorted(data_files)
    async with data_files["multiqc_general_stats.txt"].open("rb") as fh:
        general = bytes(await fh.read()).decode()
    for s in samples:
        assert s.id in general, f"{s.id} missing from MultiQC general stats"
    assert any("strand_check" in name for name in data_files), sorted(data_files)
    # Every tool's outputs were found and parsed (MultiQC writes one table per module).
    for module in ("fastqc", "cutadapt", "bbmap", "star", "samtools", "salmon", "qualimap", "rseqc", "dupradar"):
        assert any(name.lower().startswith(f"multiqc_{module}") or module in name.lower() for name in data_files), (
            f"MultiQC parsed no {module} outputs: {sorted(data_files)}"
        )

    # Published in upstream's --outdir layout.
    assert result.outdir is not None
    published = {f.path[len(result.outdir.path.rstrip("/")) :].lstrip("/") async for f in result.outdir.walk()}
    for path in [
        "fastqc/raw/WT_REP1_raw_1_fastqc.html",
        "fastqc/trim/WT_REP1_trimmed_1_val_1_fastqc.html",
        "trimgalore/WT_REP1_trimmed_1.fastq.gz_trimming_report.txt",
        "fq_lint/raw/WT_REP1.fq_lint.txt",
        "bbsplit/WT_REP1.stats.txt",
        "star_salmon/log/WT_REP1.Log.final.out",
        "star_salmon/WT_REP1.markdup.sorted.bam",
        "star_salmon/WT_REP1.markdup.sorted.bam.bai",
        "star_salmon/samtools_stats/WT_REP1.markdup.sorted.bam.flagstat",
        "star_salmon/WT_REP1/quant.sf",
        "star_salmon/salmon.merged.gene_counts.tsv",
        "star_salmon/qualimap/WT_REP1/qualimapReport.html",
        "star_salmon/rseqc/infer_experiment/WT_REP1.infer_experiment.txt",
        "star_salmon/stringtie/WT_REP1.transcripts.gtf",
        "star_salmon/bigwig/WT_REP1.forward.bigWig",
        "multiqc/star_salmon/multiqc_report.html",
    ]:
        assert path in published, f"not published: {path}"

    # Every sample (and nothing else) is a column of the merged gene matrix.
    async with m.tximport.counts_gene.open("rb") as fh:
        header = bytes(await fh.read()).decode().splitlines()[0].split("\t")
    assert sorted(header[2:]) == sorted(s.id for s in samples), header


@env.task
async def test_rnaseq_umi() -> None:
    # upstream case: tests/umi.nf.test "--umi_dedup_tool 'umitools'"
    (samplesheet, samples), fasta, gtf, transcript_fasta, additional_fasta = await asyncio.gather(
        stage_samplesheet(SAMPLESHEET),
        fixture(FASTA),
        fixture(GTF),
        fixture(TRANSCRIPT_FASTA),
        fixture(ADDITIONAL_FASTA),
    )
    salmon_index = await fixture_dir(SALMON_INDEX, "salmon")

    result = await rnaseq(
        samplesheet,
        fasta=fasta,
        gtf=gtf,
        transcript_fasta=transcript_fasta,
        additional_fasta=additional_fasta,
        salmon_index_dir=salmon_index,
        with_umi=True,
        umitools_extract_method="regex",
        umitools_bc_pattern="^(?P<umi_1>CGA.{8}){s<=2}.*",
        umitools_dedup_stats=True,
        skip_bbsplit=True,
        skip_stringtie=True,
        skip_bigwig=True,
        outdir=flyte.ctx().raw_data_path.get_random_remote_path("rnaseq_umi_outdir"),
    )

    # As upstream: every sample keeps enough reads after extraction and trimming.
    assert result.failed_trimming == {}, result.failed_trimming
    assert [r.sample for r in result.samples] == [s.id for s in samples]
    paired = {s.id for s in samples if not s.single_end}
    for r in result.samples:
        ext = r.preprocessing.umi_extract
        assert ext is not None, f"{r.sample}: no UMI extraction"
        await assert_nonempty(ext.log, label=f"{r.sample} umi_extract.log")
        assert r.markduplicates is None, f"{r.sample}: MarkDuplicates ran despite UMI deduplication"
        umi = r.umi_dedup
        assert umi is not None, f"{r.sample}: no UMI deduplication"
        assert umi.bam.path.endswith(f"{r.sample}.umi_dedup.sorted.bam"), umi.bam.path
        for label, f in [
            ("bai", umi.bai),
            ("stats", umi.stats),
            ("flagstat", umi.flagstat),
            ("idxstats", umi.idxstats),
            ("genome dedup log", umi.genome.log),
            ("transcriptome dedup log", umi.transcriptome.log),
            ("transcriptome bam", umi.transcriptome_bam),
        ]:
            await assert_nonempty(f, label=f"{r.sample} umi {label}")
        # umitools_dedup_stats: the UMI TSVs, genome and transcriptome alike.
        for side in (umi.genome, umi.transcriptome):
            for f in (side.edit_distance, side.per_umi, side.per_umi_per_position):
                assert f is not None, f"{r.sample}: missing UMI stats"
                await assert_nonempty(f, label=f"{r.sample} {f.path.rsplit('/', 1)[-1]}")
        # Upstream runs prepare-for-rsem on the paired-end samples only.
        assert (umi.prepare_for_rsem_log is not None) == (r.sample in paired), r.sample
        quant = await r.salmon.get_file("quant.sf")
        assert quant is not None, f"{r.sample}: salmon results have no quant.sf"
        await assert_nonempty(quant, label=f"{r.sample} quant.sf")

    assert result.outdir is not None
    published = {f.path[len(result.outdir.path.rstrip("/")) :].lstrip("/") async for f in result.outdir.walk()}
    for path in [
        "umitools/WT_REP1.umi_extract.log",
        "star_salmon/WT_REP1.umi_dedup.sorted.bam",
        "star_salmon/umitools/genomic_dedup_log/WT_REP1.umi_dedup.sorted.log",
        "star_salmon/umitools/prepare_for_quantification_log/WT_REP1.umi_dedup.transcriptome.filtered.prepare_for_rsem.log",
        "star_salmon/WT_REP1/quant.sf",
    ]:
        assert path in published, f"not published: {path}"
    assert not any("markdup" in p for p in published), "MarkDuplicates outputs published despite UMI deduplication"

    assert result.multiqc is not None
    data_files = [f.path.rsplit("/", 1)[-1] async for f in result.multiqc.results.walk() if "_data/" in f.path]
    assert any("umitools" in name.lower() for name in data_files), f"MultiQC parsed no UMI-tools logs: {data_files}"

    async with result.salmon.tximport.counts_gene.open("rb") as fh:
        header = bytes(await fh.read()).decode().splitlines()[0].split("\t")
    assert sorted(header[2:]) == sorted(s.id for s in samples), header


tests = [test_rnaseq_star_salmon, test_rnaseq_umi]
