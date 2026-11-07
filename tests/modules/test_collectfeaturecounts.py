"""Tests for flyte_bio.modules.collectfeaturecounts.

The cases and expected values are upstream's custom/collectfeaturecounts
tests (modules commit 9a1cc9f). As upstream, the counts table is checked line
by line after decompression, and the versions file by its md5. Upstream's
inline inputs are written by the tests; its failure case checks that the task
fails with the R script's message.
"""

import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.collectfeaturecounts import collect_featurecounts
from tests.framework import assert_md5, env, fixture

VERSIONS_MD5 = "c7e1062ade0a362673244b19bcf210cf"
HEADER = "# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t"
COLUMNS = "orf\tchr\tstart\tend\tstrand\tlength\tsample\tcount\ttpm"


async def write_table(name: str, content: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(content)
    return await File.from_local(str(path))


async def lines_gzip(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode().splitlines()


@env.task
async def test_collectfeaturecounts_two_samples() -> None:
    # upstream case: custom/collectfeaturecounts "two samples"
    tables = [
        await fixture("generic/tsv/sample1.subread.featureCounts.tsv"),
        await fixture("generic/tsv/sample2.subread.featureCounts.tsv"),
    ]
    result = await collect_featurecounts(tables, prefix="test")
    assert result.counts.name == "test.counts.tsv.gz", f"wrote {result.counts.name}"
    expected = [
        COLUMNS,
        "orf1ab\tMT192765.1;MT192765.1\t259;13461\t13461;21545\t+;+\t21287\tSAMPLE1\t38\t204755.733883",
        "S\tMT192765.1\t21556\t25374\t+\t3819\tSAMPLE1\t7\t210239.985324",
        "ORF3a\tMT192765.1\t25386\t26210\t+\t825\tSAMPLE1\t1\t139031.429256",
        "M\tMT192765.1\t26516\t27181\t+\t666\tSAMPLE1\t1\t172223.617322",
        "N\tMT192765.1\t28267\t29523\t+\t1257\tSAMPLE1\t3\t273749.234215",
        "orf1ab\tMT192765.1;MT192765.1\t259;13461\t13461;21545\t+;+\t21287\tSAMPLE2\t52\t125022.072833",
        "S\tMT192765.1\t21556\t25374\t+\t3819\tSAMPLE2\t4\t53605.351066",
        "E\tMT192765.1\t26238\t26462\t+\t225\tSAMPLE2\t2\t454930.746049",
        "N\tMT192765.1\t28267\t29523\t+\t1257\tSAMPLE2\t9\t366441.830052",
    ]
    lines = await lines_gzip(result.counts)
    assert lines == expected, f"counts: {lines}"
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


@env.task
async def test_collectfeaturecounts_sorted_bam_in_sample_name() -> None:
    # upstream case: custom/collectfeaturecounts "sample name containing sorted and bam"
    table = await write_table(
        "rawsorted_bam.featureCounts.tsv",
        "# Program:featureCounts\n"
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\trawsorted_bam.sorted.bam\n"
        "gene1\tcontig1\t1\t300\t+\t300\t10\n",
    )
    result = await collect_featurecounts([table], prefix="test")
    sample = (await lines_gzip(result.counts))[1].split("\t")[6]
    assert sample == "rawsorted_bam", f"sample: {sample}"


@env.task
async def test_collectfeaturecounts_first_sample_without_counts() -> None:
    # upstream case: custom/collectfeaturecounts "first sample without counts"
    tables = [
        await write_table(
            "a.featureCounts.tsv", HEADER + "A.sorted.bam\ng1\tc\t1\t100\t+\t100\t0\ng2\tc\t200\t399\t+\t200\t0\n"
        ),
        await write_table(
            "b.featureCounts.tsv", HEADER + "B.sorted.bam\ng1\tc\t1\t100\t+\t100\t3\ng2\tc\t200\t399\t+\t200\t0\n"
        ),
    ]
    result = await collect_featurecounts(tables, prefix="test")
    lines = await lines_gzip(result.counts)
    assert lines == [COLUMNS, "g1\tc\t1\t100\t+\t100\tB\t3\t1000000"], f"counts: {lines}"
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


@env.task
async def test_collectfeaturecounts_same_sample_in_two_files_fails() -> None:
    # upstream case: custom/collectfeaturecounts "the same sample in two files fails"
    tables = [
        await write_table("a.featureCounts.tsv", HEADER + "A.sorted.bam\ng1\tc\t1\t100\t+\t100\t3\n"),
        await write_table("b.featureCounts.tsv", HEADER + "A.sorted.bam\ng2\tc\t200\t399\t+\t200\t4\n"),
    ]
    try:
        await collect_featurecounts(tables, prefix="test")
    except Exception as e:
        assert "is in more than one input file" in str(e), f"the error lacks the duplicate-sample message: {e}"
        return
    raise AssertionError("collect_featurecounts accepted a sample in two files")


tests = [
    test_collectfeaturecounts_two_samples,
    test_collectfeaturecounts_sorted_bam_in_sample_name,
    test_collectfeaturecounts_first_sample_without_counts,
    test_collectfeaturecounts_same_sample_in_two_files_fails,
]
