"""Tests for flyte_bio.modules.rustqc.

The case is upstream's rustqc test (modules commit 9a1cc9f, RustQC 0.2.1),
with upstream's test ``ext.args`` (``--skip-dup-check --stranded forward``)
and its assertions: md5s of the reproducible tables (featureCounts, preseq,
and RSeQC's infer_experiment, TIN and summaries), the file names of the other
outputs, and upstream's content checks on samtools flagstat and idxstats,
Qualimap and dupRadar.

The case runs on 2 CPUs, as upstream's did: featureCounts.tsv's header records
the command line, thread count (``-t``) included.
"""

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.rustqc import rustqc_cmd
from tests.framework import assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"
FEATURECOUNTS_MD5S = {
    "test.biotype_counts.tsv": "68744ae62d861b2f5fbe83b2d1e1ac92",
    "test.biotype_counts_mqc.tsv": "0bb6b11917b544c9690e0c9104f044d6",
    "test.biotype_counts_rrna_mqc.tsv": "31c68d05ffe90d1e3eb2cb2ef9a7b4d9",
    "test.featureCounts.biotype.tsv.summary": "45ccc3a9a51380ab0607d1e87bbd9db6",
    "test.featureCounts.tsv": "7613a7c0da252ea99f69ac7c33cd8e7c",
    "test.featureCounts.tsv.summary": "45ccc3a9a51380ab0607d1e87bbd9db6",
}
PRESEQ_MD5S = {"test.lc_extrap.txt": "deefbb79ddd19e5d05f4f757ed2f2f13"}
RSEQC_MD5S = {
    "test.infer_experiment.txt": "58428af8e55b2a464a247650ca926908",
    "test.inner_distance_summary.txt": "d8e829669dca96f3fe485f026258726e",
    "test.junctionSaturation_summary.txt": "622a2ebe648988edddb9aab2f6f11678",
    "test.summary.txt": "0c1c79cc34b229369eb027c7bafada3b",
    "test.tin.xls": "e5699ed57f9b0513f9e7115de2dd855a",
}
DUPRADAR_NAMES = [
    "test_dupMatrix.txt",
    "test_duprateExpBoxplot.png",
    "test_duprateExpBoxplot.svg",
    "test_expressionHist.png",
    "test_expressionHist.svg",
]
SAMTOOLS_NAMES = ["test.flagstat", "test.idxstats", "test.stats"]
QUALIMAP_NAMES = [
    "Coverage Profile Along Genes (High).png",
    "Coverage Profile Along Genes (High).svg",
    "Coverage Profile Along Genes (Low).png",
    "Coverage Profile Along Genes (Low).svg",
    "Coverage Profile Along Genes (Total).png",
    "Coverage Profile Along Genes (Total).svg",
    "Reads Genomic Origin.png",
    "Reads Genomic Origin.svg",
    "Transcript coverage histogram.png",
    "Transcript coverage histogram.svg",
    "coverage_profile_along_genes_(high).txt",
    "coverage_profile_along_genes_(low).txt",
    "coverage_profile_along_genes_(total).txt",
    "qualimapReport.html",
    "rnaseq_qc_results.txt",
]
RSEQC_OTHER_NAMES = [
    "test.DupRate_plot.png",
    "test.DupRate_plot.r",
    "test.DupRate_plot.svg",
    "test.bam_stat.txt",
    "test.inner_distance.txt",
    "test.inner_distance_freq.txt",
    "test.inner_distance_mean.txt",
    "test.inner_distance_plot.png",
    "test.inner_distance_plot.r",
    "test.inner_distance_plot.svg",
    "test.junction.Interact.bed",
    "test.junction.bed",
    "test.junction.xls",
    "test.junctionSaturation_plot.png",
    "test.junctionSaturation_plot.r",
    "test.junctionSaturation_plot.svg",
    "test.junction_annotation.log",
    "test.junction_plot.r",
    "test.pos.DupRate.xls",
    "test.read_distribution.txt",
    "test.seq.DupRate.xls",
    "test.splice_events.png",
    "test.splice_events.svg",
    "test.splice_junction.png",
    "test.splice_junction.svg",
]


async def rustqc_files(results: Dir, part: str) -> dict[str, File]:
    """Every file under ``<results>/<part>/``, recursively, by name."""
    files: dict[str, File] = {}
    prefix = results.path.rstrip("/") + f"/{part}/"
    async for f in results.walk():
        if f.path.startswith(prefix):
            name = f.path.rsplit("/", 1)[-1]
            assert name not in files, f"{part} holds two files named {name}"
            files[name] = f
    return files


async def rustqc_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_rustqc_paired_end_bam() -> None:
    # upstream case: rustqc "homo_sapiens paired-end [bam]"
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    task = rustqc_cmd.override(
        resources=flyte.Resources(cpu=2, memory="16Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=2")
    )
    output = await task(
        bam=await fixture(BAM),
        bai=await fixture(BAM + ".bai"),
        gtf=await fixture("genomics/homo_sapiens/genome/genome.gtf"),
        paired=True,
        prefix="test",
        args="--skip-dup-check --stranded forward",
    )
    results = Dir.from_existing_remote(output.path.rstrip("/") + "/test")
    featurecounts = await rustqc_files(results, "featurecounts")
    preseq = await rustqc_files(results, "preseq")
    rseqc = await rustqc_files(results, "rseqc")
    dupradar = await rustqc_files(results, "dupradar")
    samtools = await rustqc_files(results, "samtools")
    qualimap = await rustqc_files(results, "qualimap")
    for part, files, md5s in [("featurecounts", featurecounts, FEATURECOUNTS_MD5S), ("preseq", preseq, PRESEQ_MD5S)]:
        assert sorted(files) == sorted(md5s), f"{part} holds {sorted(files)}"
        for name, md5 in md5s.items():
            await assert_md5(files[name], md5, label=name)
    for name, md5 in RSEQC_MD5S.items():
        assert name in rseqc, f"rseqc lacks {name}"
        await assert_md5(rseqc[name], md5, label=name)
    for part, names, expected in [
        ("dupradar", sorted(dupradar), DUPRADAR_NAMES),
        ("samtools", sorted(samtools), SAMTOOLS_NAMES),
        ("qualimap", sorted(qualimap), QUALIMAP_NAMES),
        ("rseqc", sorted(set(rseqc) - set(RSEQC_MD5S)), RSEQC_OTHER_NAMES),
    ]:
        assert names == expected, f"{part} holds {names}, expected {expected}"
    assert "5644 + 0 in total" in await rustqc_text(samtools["test.flagstat"]), "flagstat lacks 5644 reads"
    assert "chr22\t40001\t5642" in await rustqc_text(samtools["test.idxstats"]), "idxstats lacks 5642 reads on chr22"
    qualimap_results = await rustqc_text(qualimap["rnaseq_qc_results.txt"])
    assert "total alignments = 5,642" in qualimap_results, "qualimap lacks 5,642 alignments"
    assert "exonic =  1,328 (23.55%)" in qualimap_results, "qualimap lacks 1,328 (23.55%) exonic"
    assert "ENSG00000239435\t3389\t664" in await rustqc_text(dupradar["test_dupMatrix.txt"]), (
        "dupMatrix lacks 664 reads"
    )


tests = [test_rustqc_paired_end_bam]
