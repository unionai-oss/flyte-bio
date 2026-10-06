"""Tests for flyte_bio.modules.qualimap.

Mirrors upstream qualimap/rnaseq: the report HTML exists and
``rnaseq_qc_results.txt`` matches the upstream snapshot md5 (the file
records the input names, so fixtures keep their original basenames).
"""

from flyte_bio.modules.qualimap import qualimap_rnaseq
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


tests = [test_qualimap_rnaseq]
