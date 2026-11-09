"""Tests for flyte_bio.modules.qualimap.

Mirrors upstream qualimap/rnaseq: the report HTML exists and
``rnaseq_qc_results.txt`` matches the upstream snapshot md5 (the file
records the input names, so fixtures keep their original basenames).

The bamqc case is upstream's qualimap/bamqc test (modules commit 9a1cc9f,
Qualimap 2.3), with its test config's ``-p non-strand-specific``. Upstream
snapshots only the report's file names (its contents record dates and
paths); this test checks the same names and that ``genome_results.txt``
summarises the staged BAM.

The bamqccram case is upstream's qualimap/bamqccram test (same commit), with
its test config's ``-p non-strand-specific``: the same report files, from a
CRAM streamed to Qualimap as ``/dev/stdin``.
"""

from flyte.io import Dir

from flyte_bio.modules.qualimap import qualimap_bamqc, qualimap_bamqccram, qualimap_rnaseq
from tests.framework import assert_md5, env, fixture


@env.task
async def test_qualimap_rnaseq() -> None:
    # upstream case: qualimap/rnaseq "homo_sapiens [bam]"
    bam = await fixture("genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam")
    gtf = await fixture("genomics/homo_sapiens/genome/genome.gtf")
    results = await qualimap_rnaseq(bam, gtf, prefix="test", single_end=False)
    assert await results.get_file("qualimapReport.html") is not None, "no qualimapReport.html"
    report = await results.get_file("rnaseq_qc_results.txt")
    assert report is not None, "no rnaseq_qc_results.txt"
    await assert_md5(report, "b77878cac45beaa79a892af54aad2da3", label="rnaseq_qc_results.txt")


BAMQC_CSS = [
    "agogo.css",
    "ajax-loader.gif",
    "basic.css",
    "bgfooter.png",
    "bgtop.png",
    "comment-bright.png",
    "comment-close.png",
    "comment.png",
    "doctools.js",
    "down-pressed.png",
    "down.png",
    "file.png",
    "jquery.js",
    "minus.png",
    "plus.png",
    "pygments.css",
    "qualimap_logo_small.png",
    "report.css",
    "searchtools.js",
    "underscore.js",
    "up-pressed.png",
    "up.png",
    "websupport.js",
]
BAMQC_IMAGES = [
    "genome_coverage_0to50_histogram.png",
    "genome_coverage_across_reference.png",
    "genome_coverage_histogram.png",
    "genome_coverage_quotes.png",
    "genome_gc_content_per_window.png",
    "genome_homopolymer_indels.png",
    "genome_insert_size_across_reference.png",
    "genome_insert_size_histogram.png",
    "genome_mapping_quality_across_reference.png",
    "genome_mapping_quality_histogram.png",
    "genome_reads_clipping_profile.png",
    "genome_reads_content_per_read_position.png",
    "genome_uniq_read_starts_histogram.png",
]
BAMQC_RAW_DATA = [
    "coverage_across_reference.txt",
    "coverage_histogram.txt",
    "duplication_rate_histogram.txt",
    "genome_fraction_coverage.txt",
    "homopolymer_indels.txt",
    "insert_size_across_reference.txt",
    "insert_size_histogram.txt",
    "mapped_reads_clipping_profile.txt",
    "mapped_reads_gc-content_distribution.txt",
    "mapped_reads_nucleotide_content.txt",
    "mapping_quality_across_reference.txt",
    "mapping_quality_histogram.txt",
]


async def check_bamqc_report(results: Dir, bam_name: str) -> None:
    root = results.path.rstrip("/") + "/"
    found = sorted([f.path.removeprefix(root) async for f in results.walk()])
    expected = sorted(
        [f"css/{n}" for n in BAMQC_CSS]
        + ["genome_results.txt", "qualimapReport.html"]
        + [f"images_qualimapReport/{n}" for n in BAMQC_IMAGES]
        + [f"raw_data_qualimapReport/{n}" for n in BAMQC_RAW_DATA]
    )
    assert found == expected, (
        f"bamqc report files differ: missing {set(expected) - set(found)}, extra {set(found) - set(expected)}"
    )
    genome_results = await results.get_file("genome_results.txt")
    assert genome_results is not None, "no genome_results.txt"
    async with genome_results.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    for wanted in [
        f"bam file = {bam_name}",
        ">>>>>>> Reference",
        "number of reads",
        "mean coverageData",
    ]:
        assert wanted in text, f"genome_results.txt lacks {wanted!r}"


@env.task
async def test_qualimap_bamqc() -> None:
    # upstream case: qualimap/bamqc "homo_sapiens [bam]"
    bam = await fixture("genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam")
    results = await qualimap_bamqc(bam, prefix="test", single_end=False, args="-p non-strand-specific")
    await check_bamqc_report(results, "test.paired_end.sorted.bam")


@env.task
async def test_qualimap_bamqccram() -> None:
    # upstream case: qualimap/bamqccram "homo_sapiens [cram]"
    cram = "genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram"
    genome = "genomics/homo_sapiens/genome/genome.fasta"
    results = await qualimap_bamqccram(
        await fixture(cram),
        await fixture(cram + ".crai"),
        await fixture(genome),
        await fixture(genome + ".fai"),
        prefix="test",
        single_end=False,
        args="-p non-strand-specific",
    )
    await check_bamqc_report(results, "/dev/stdin")


tests = [test_qualimap_rnaseq, test_qualimap_bamqc, test_qualimap_bamqccram]
