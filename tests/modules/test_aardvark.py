"""Tests for flyte_bio.modules.aardvark.

The cases are upstream's aardvark/compare tests (modules commit 9a1cc9f,
Aardvark 0.10.5), plain and stratified by ``multi_intervals.bed.gz``: nft-vcf's
summary of the labelled truth and query VCFs (contigs, samples, records) and
their names, and the md5s of the region tables (decompressed, as nf-test
hashes them), the summary and the settings JSON. The settings record the
thread count, so the tests run on the 4 CPUs upstream's snapshot ran on.
"""

import asyncio
import gzip
import tempfile
from pathlib import Path

import flyte
from flyte.io import File

from flyte_bio.modules.aardvark import AardvarkCompareResult, aardvark_compare_cmd, aardvark_compare_result
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
QUERY = HS + "illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz"
TRUTH = HS + "illumina/gatk/haplotypecaller_calls/test2_haplotc.vcf.gz"
CHR21 = HS + "genome/chr21/sequence/"


async def aardvark_vcf_summary(vcf: File) -> tuple[list[str], int, int]:
    """A gzipped VCF's contigs (in order of first record), sample count and record count."""
    async with vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    samples = len(next(x for x in lines if x.startswith("#CHROM")).split("\t")) - 9
    records = [x.split("\t", 1)[0] for x in lines if not x.startswith("#")]
    return list(dict.fromkeys(records)), samples, len(records)


async def aardvark_case(stratified: bool) -> AardvarkCompareResult:
    query, query_tbi, truth, truth_tbi, regions, fasta = await asyncio.gather(
        fixture(QUERY),
        fixture(QUERY + ".tbi"),
        fixture(TRUTH),
        fixture(TRUTH + ".tbi"),
        fixture(CHR21 + "multi_intervals.bed"),
        fixture(CHR21 + "genome.fasta"),
    )
    tsv, beds = None, None
    if stratified:
        path = Path(tempfile.mkdtemp()) / "stratification.tsv"
        path.write_text("multi_intervals\tmulti_intervals.bed.gz\n")
        tsv, beds = await File.from_local(str(path)), [await fixture(CHR21 + "multi_intervals.bed.gz")]
    # The settings record the thread count; upstream's snapshot ran on 4 CPUs. The cache key ignores resources, so
    # salt it to keep these runs apart from default ones.
    task = aardvark_compare_cmd.override(
        resources=flyte.Resources(cpu=4, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=4")
    )
    results = await task(
        query_vcf=query,
        query_tbi=query_tbi,
        truth_vcf=truth,
        truth_tbi=truth_tbi,
        regions_bed=regions,
        fasta=fasta,
        stratification_tsv=tsv,
        stratification_beds=beds or [],
        prefix="test",
        args="",
    )
    return await aardvark_compare_result(results, "test")


async def check(result: AardvarkCompareResult, summary_md5: str, runinfo_md5: str) -> None:
    assert result.labelled_truth.name == "test.truth.vcf.gz", result.labelled_truth.name
    assert result.labelled_query.name == "test.query.vcf.gz", result.labelled_query.name
    truth = await aardvark_vcf_summary(result.labelled_truth)
    assert truth == (["chr21"], 1, 737), f"labelled truth (contigs, samples, records) {truth}"
    query = await aardvark_vcf_summary(result.labelled_query)
    assert query == ([], 1, 0), f"labelled query (contigs, samples, records) {query}"
    await assert_gunzipped_md5(result.region_sequences, "9252f664552feb5e6a406f3757513411", label="region sequences")
    await assert_gunzipped_md5(result.region_summary, "661d4e2335bdf0071712077ac21ee4f8", label="region summary")
    await assert_md5(result.summary, summary_md5, label="summary")
    await assert_md5(result.runinfo, runinfo_md5, label="cli settings")


@env.task
async def test_aardvark_compare_chr21() -> None:
    await check(await aardvark_case(False), "02c4580cae7bbbcbc6f08f6c3c33f229", "3bb5f03eb62d80bb878a85df7f89f89b")


@env.task
async def test_aardvark_compare_stratification() -> None:
    await check(await aardvark_case(True), "dae4540f9abb90673d377506ee7ecd67", "758b7872e9ba675d7b4dbe3dad9bcf6d")


tests = [test_aardvark_compare_chr21, test_aardvark_compare_stratification]
