"""Tests for flyte_bio.modules.happy.

The cases are upstream's happy/happy, happy/sompy and happy/report tests
(modules commit 9a1cc9f, hap.py 0.3.15, ga4gh-happy-report 0.1.2):

- hap.py, within regions, with a false-positive BED, stratified and plain:
  the md5s of the summary, extended and ROC tables (decompressed, as nf-test
  hashes them), the names of the metrics, run-info and index files, and the
  annotated VCF's ``##fileformat=VCFv4.1`` line;
- som.py with the test config's ``--feature-table hcc.mutect.snv``, without
  and with a BAM: the stats table's md5 and the metrics' and features' names;
- the report of one ROC table and of two labelled ones: the HTML's name.
"""

import asyncio
import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.happy import HappyResult, happy_happy, happy_report, happy_sompy
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
HAPLOTC = HS + "illumina/gatk/haplotypecaller_calls/test2_haplotc"
MUTECT2 = HS + "illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz"
CHR21 = HS + "genome/chr21/sequence/"
SOMPY_ARGS = "--feature-table hcc.mutect.snv"


async def happy_case(prefix: str = "test", **optional: File | list[File] | None) -> HappyResult:
    query, truth, fasta, fai = await asyncio.gather(
        fixture(HAPLOTC + ".vcf.gz"),
        fixture(HAPLOTC + ".ann.vcf.gz"),
        fixture(CHR21 + "genome.fasta"),
        fixture(CHR21 + "genome.fasta.fai"),
    )
    return await happy_happy(query, truth, fasta, fai, prefix=prefix, **optional)


async def check_happy(result: HappyResult, md5s: dict[str, str]) -> None:
    assert result.vcf is not None and result.tbi is not None, "hap.py wrote no VCF"
    async with result.vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    assert "##fileformat=VCFv4.1" in lines, "no ##fileformat=VCFv4.1 line"
    for key, md5 in md5s.items():
        file = getattr(result, key)
        if file.name.endswith(".gz"):
            await assert_gunzipped_md5(file, md5, label=key)
        else:
            await assert_md5(file, md5, label=key)
    names = (result.metrics_json.name, result.runinfo.name, result.tbi.name)
    assert names == ("test.metrics.json.gz", "test.runinfo.json", "test.vcf.gz.tbi"), names


def happy_md5s(summary: str, extended: str, roc_all: str, indel: str, snp: str, indel_pass: str, snp_pass: str):
    return {
        "summary_csv": summary,
        "extended_csv": extended,
        "roc_all_csv": roc_all,
        "roc_indel_locations_csv": indel,
        "roc_snp_locations_csv": snp,
        "roc_indel_locations_pass_csv": indel_pass,
        "roc_snp_locations_pass_csv": snp_pass,
    }


@env.task
async def test_happy_happy_bed() -> None:
    result = await happy_case(regions_bed=await fixture(CHR21 + "multi_intervals.bed"))
    md5s = happy_md5s(
        "9365186b70c0bfb71630ba16d7b5089b",
        "494f93d3a2f5ea5af9debece4b0e4055",
        "e293f78fb6b7a1ef5fdfa95c2fbc371a",
        "e7035c81472ea713ebdaa4c7f6fbc43b",
        "c3e63e548b9937192b19841a10f1f199",
        "c7b59f00a40bc2ce1b8d8e0c04caa7a8",
        "a2f8cccb01841be15163671b1d356962",
    )
    await check_happy(result, md5s)


@env.task
async def test_happy_happy_false_positives() -> None:
    result = await happy_case(false_positives_bed=await fixture(CHR21 + "multi_intervals.bed"))
    md5s = happy_md5s(
        "8a302b92d08de8d68509e3adc4bb9d4d",
        "3435c0b5d2adec043e72de9cf86deba0",
        "6d05ef088fb31469f17c379603c00a25",
        "4e670c9fd9e27ef3b2fff7cb1064ef22",
        "c6d6e8b69864927767b755309eca3b58",
        "b2e30f99571d167a938ce0c35d6e9202",
        "1a74e0cfe988b7278adb49c1058d1b9b",
    )
    await check_happy(result, md5s)


@env.task
async def test_happy_happy_stratification() -> None:
    path = Path(tempfile.mkdtemp()) / "stratification.tsv"
    path.write_text("multi_intervals\tmulti_intervals.bed.gz\n")
    result = await happy_case(
        stratification_tsv=await File.from_local(str(path)),
        stratification_beds=[await fixture(CHR21 + "multi_intervals.bed.gz")],
    )
    md5s = happy_md5s(
        "d802ef043c6b504e7a103e8dcd4ba300",
        "5ffe4d6bf90d05bd65024265adc3938c",
        "2736ef78791a0aa8a31f611ed0e67cb4",
        "e7035c81472ea713ebdaa4c7f6fbc43b",
        "12270ccd8b277097ac899b602de1fea9",
        "c7b59f00a40bc2ce1b8d8e0c04caa7a8",
        "16466cda713d728d91e801392b514158",
    )
    await check_happy(result, md5s)


@env.task
async def test_happy_happy_nobed() -> None:
    md5s = happy_md5s(
        "d802ef043c6b504e7a103e8dcd4ba300",
        "c7d1851ffb1aef2b2e3ca4aaf94dade3",
        "f1ccbf86653374a9c09a20f5dec3bc50",
        "e7035c81472ea713ebdaa4c7f6fbc43b",
        "12270ccd8b277097ac899b602de1fea9",
        "c7b59f00a40bc2ce1b8d8e0c04caa7a8",
        "16466cda713d728d91e801392b514158",
    )
    await check_happy(await happy_case(), md5s)


async def sompy_case(bams: list[File], stats_md5: str) -> None:
    query, truth, regions, fasta, fai = await asyncio.gather(
        fixture(MUTECT2),
        fixture(HAPLOTC + ".vcf.gz"),
        fixture(CHR21 + "multi_intervals.bed"),
        fixture(CHR21 + "genome.fasta"),
        fixture(CHR21 + "genome.fasta.fai"),
    )
    result = await happy_sompy(query, truth, fasta, fai, prefix="test", regions_bed=regions, bams=bams, args=SOMPY_ARGS)
    await assert_md5(result.stats, stats_md5, label="stats")
    assert result.features is not None, "som.py wrote no feature table"
    names = (result.metrics.name, result.features.name)
    assert names == ("test.metrics.json", "test.features.csv"), names


@env.task
async def test_happy_sompy_vcf() -> None:
    await sompy_case([], "8b35dd00bf7b5bd697a05d73fb8c0816")


@env.task
async def test_happy_sompy_vcf_bam() -> None:
    bam = await fixture(HS + "illumina/bam/test2.paired_end.sorted.bam")
    await sompy_case([bam], "5e9aae3e92641f1934c10fe88f250e59")


@env.task
async def test_happy_report_single_roc() -> None:
    roc = (await happy_case()).roc_all_csv
    report = await happy_report([roc], ["gatk4_vcfeval"], prefix="test")
    assert report.name == "test.html", report.name


@env.task
async def test_happy_report_multiple_labelled_roc() -> None:
    plain, false_positives = await asyncio.gather(
        happy_case(), happy_case(prefix="test_fp", false_positives_bed=await fixture(CHR21 + "multi_intervals.bed"))
    )
    report = await happy_report(
        [plain.roc_all_csv, false_positives.roc_all_csv], ["gatk4_vcfeval", "gatk4_happy"], prefix="compare"
    )
    assert report.name == "compare.html", report.name


tests = [
    test_happy_happy_bed,
    test_happy_happy_false_positives,
    test_happy_happy_stratification,
    test_happy_happy_nobed,
    test_happy_sompy_vcf,
    test_happy_sompy_vcf_bam,
    test_happy_report_single_roc,
    test_happy_report_multiple_labelled_roc,
]
