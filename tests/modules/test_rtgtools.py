"""Tests for flyte_bio.modules.rtgtools.

The cases are upstream's rtgtools/format tests (modules commit 9a1cc9f, RTG
Tools 3.13), whose snapshots list the files of the SDF directory: sarscov2's
genome FASTA, its paired FASTQs (``--format fastq``) and its single-end BAM
with a read group (``--format sam-se``).

As upstream, the FASTQ case passes both files but leaves ``single_end`` at its
default, so only ``test_1.fastq.gz`` is formatted.

The rtgtools/vcfeval, rtgtools/bndeval, rtgtools/cnveval and
rtgtools/svdecompose cases are upstream's at the same commit, with their
snapshot md5s (nf-test hashes ``.gz`` outputs decompressed, indexes as they
are). vcfeval's reference is upstream's chr21 SDF fixture, linked under its
upstream name, ``test`` (upstream's untar names its output after the test's
id), which vcfeval records; vcfeval and cnveval run
on 2 CPUs, as upstream's snapshots did, since they record ``--threads``; bndeval's query is
svdecompose's output, as upstream.
"""

import asyncio
import tempfile
from pathlib import Path

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.rtgtools import (
    cnveval_result,
    rtgtools_bndeval,
    rtgtools_cnveval_cmd,
    rtgtools_format,
    rtgtools_svdecompose,
    rtgtools_vcfeval_cmd,
    vcfeval_result,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture, fixture_dir

HS = "genomics/homo_sapiens/"
HAPLOTC = HS + "illumina/gatk/haplotypecaller_calls/test2_haplotc"
SIMULATED_SV = HS + "illumina/vcf/chr21/simulated_sv"
EMPTY_TBI = "092a7a3162e7cff25d273525751eb284"
# vcfeval and cnveval record --threads in their headers; upstream's snapshots ran on 2 CPUs. The cache key ignores
# resources, so salt it to keep these runs apart from default ones.
TWO_CPUS = {"resources": flyte.Resources(cpu=2, memory="8Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=2")}

COMMON = [
    "done",
    "format.log",
    "mainIndex",
    "nameIndex0",
    "namedata0",
    "namepointer0",
    "progress",
    "seqdata0",
    "seqpointer0",
    "sequenceIndex0",
    "summary.txt",
]
SUFFIX = ["suffixIndex0", "suffixdata0", "suffixpointer0"]


async def sdf_files(sdf: Dir) -> list[str]:
    return sorted(f.name for f in await sdf.list_files())


async def check_sdf(sdf: Dir, expected: list[str]) -> None:
    assert sdf.path.rstrip("/").endswith("/test.sdf"), sdf.path
    names = await sdf_files(sdf)
    assert names == sorted(expected), f"SDF files {names}"


@env.task
async def test_rtgtools_format_fasta() -> None:
    sdf = await rtgtools_format(await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test")
    await check_sdf(sdf, COMMON + SUFFIX)


@env.task
async def test_rtgtools_format_fastqs() -> None:
    sdf = await rtgtools_format(
        await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz"),
        await fixture("genomics/sarscov2/illumina/fastq/test_2.fastq.gz"),
        prefix="test",
        args="--format fastq",
    )
    await check_sdf(sdf, [*COMMON, "qualitydata0", *SUFFIX])


@env.task
async def test_rtgtools_format_bam_rg() -> None:
    rg = Path(tempfile.mkdtemp()) / "sam_rg.txt"
    rg.write_text("@RG\tID:READGROUP1\tSM:SAMPLE\tPL:ILLUMINA")
    sdf = await rtgtools_format(
        await fixture("genomics/sarscov2/illumina/bam/test.single_end.bam"),
        sam_rg=await File.from_local(str(rg)),
        prefix="test",
        args="--format sam-se",
    )
    await check_sdf(sdf, [*COMMON, "qualitydata0"])


async def check_md5s(files: dict[str, tuple[File, str]]) -> None:
    """Each named file's md5 (decompressed for ``.gz``), as nf-test snapshots it."""
    for label, (file, md5) in files.items():
        if file.name.endswith(".gz"):
            await assert_gunzipped_md5(file, md5, label=label)
        else:
            await assert_md5(file, md5, label=label)


async def vcfeval_sdf() -> Dir:
    return await fixture_dir(HS + "genome/chr21/sequence/genome_sdf.tar.gz", "genome_sdf")


@env.task
async def test_rtgtools_vcfeval_all_inputs() -> None:
    query, query_tbi, truth, truth_tbi, truth_bed, regions_bed, sdf = await asyncio.gather(
        fixture(HAPLOTC + ".vcf.gz"),
        fixture(HAPLOTC + ".vcf.gz.tbi"),
        fixture(HAPLOTC + ".ann.vcf.gz"),
        fixture(HAPLOTC + ".ann.vcf.gz.tbi"),
        fixture(HS + "genome/genome.bed"),
        fixture(HS + "genome/chr21/sequence/multi_intervals.bed"),
        vcfeval_sdf(),
    )
    results = await rtgtools_vcfeval_cmd.override(**TWO_CPUS)(
        query_vcf=query,
        query_tbi=query_tbi,
        truth_vcf=truth,
        truth_tbi=truth_tbi,
        truth_bed=truth_bed,
        regions_bed=regions_bed,
        sdf=sdf,
        sdf_name="test",
        prefix="test",
        args="",
    )
    result = await vcfeval_result(results, "test")
    await check_md5s(
        {
            "tp": (result.tp_vcf, "5171021307097220337dbcaccc860495"),
            "tp tbi": (result.tp_tbi, EMPTY_TBI),
            "fn": (result.fn_vcf, "fc419367818700d47df073615aeb9077"),
            "fn tbi": (result.fn_tbi, EMPTY_TBI),
            "fp": (result.fp_vcf, "5171021307097220337dbcaccc860495"),
            "fp tbi": (result.fp_tbi, EMPTY_TBI),
            "baseline": (result.baseline_vcf, "fc419367818700d47df073615aeb9077"),
            "baseline tbi": (result.baseline_tbi, EMPTY_TBI),
            "snp roc": (result.snp_roc, "6c8977c744a08c90123a21f22c0f4e54"),
            "non-snp roc": (result.non_snp_roc, "574afefc756c4102bb3aefb0cde19dff"),
            "weighted roc": (result.weighted_roc, "c0155caab6b8a314c303ca1b481ba1af"),
            "summary": (result.summary, "f4c8df93c8bdab603036bbc27b4a28c3"),
            "phasing": (result.phasing, "31988234bee208cacb3de90dabe1797f"),
        }
    )


@env.task
async def test_rtgtools_vcfeval_no_indexes() -> None:
    query, truth, sdf = await asyncio.gather(
        fixture(HAPLOTC + ".vcf.gz"), fixture(HAPLOTC + ".ann.vcf.gz"), vcfeval_sdf()
    )
    results = await rtgtools_vcfeval_cmd.override(**TWO_CPUS)(
        query_vcf=query,
        query_tbi=None,
        truth_vcf=truth,
        truth_tbi=None,
        truth_bed=None,
        regions_bed=None,
        sdf=sdf,
        sdf_name="test",
        prefix="test",
        args="",
    )
    result = await vcfeval_result(results, "test")
    await check_md5s(
        {
            "tp": (result.tp_vcf, "5125ee41457c9d93f46b19e32788edb4"),
            "tp tbi": (result.tp_tbi, "a0e9ac2d38c04bd591ab8f857c5c9133"),
            "fn": (result.fn_vcf, "df96e4e4014cdb3050cb6f221f0cdca9"),
            "fn tbi": (result.fn_tbi, EMPTY_TBI),
            "fp": (result.fp_vcf, "d4bfa2c7271351ca19589f0f57f210b7"),
            "fp tbi": (result.fp_tbi, EMPTY_TBI),
            "baseline": (result.baseline_vcf, "920af25c3c18a438b11440702562fa35"),
            "baseline tbi": (result.baseline_tbi, "95938320b425e28cf06c45ab45ad0360"),
            "snp roc": (result.snp_roc, "7906ac0bdf54d34ac9b52bfd101cddfd"),
            "non-snp roc": (result.non_snp_roc, "0917d6986ae90818f431c316b1cbea27"),
            "weighted roc": (result.weighted_roc, "599331bcba15312abc402bdffc8d2e30"),
            "summary": (result.summary, "f33feb32f84958fb931063044fba369b"),
            "phasing": (result.phasing, "133677dbd8be657439ea2b03fdfb8795"),
        }
    )


@env.task
async def test_rtgtools_svdecompose_vcf() -> None:
    vcf, tbi = await asyncio.gather(fixture(SIMULATED_SV + ".vcf.gz"), fixture(SIMULATED_SV + ".vcf.gz.tbi"))
    result = await rtgtools_svdecompose(vcf, prefix="test", tbi=tbi)
    await check_md5s(
        {
            "vcf": (result.vcf, "609e350f55181c32323cc83c6c96a3ac"),
            "tbi": (result.tbi, "7ee21126b9f620992b0bcb5efe528f70"),
        }
    )


@env.task
async def test_rtgtools_bndeval_regions() -> None:
    vcf, tbi, regions = await asyncio.gather(
        fixture(SIMULATED_SV + ".vcf.gz"), fixture(SIMULATED_SV + ".vcf.gz.tbi"), fixture(HS + "genome/genome.bed")
    )
    # As upstream's channel join: svdecompose's output is the query, the original VCF the truth.
    decomposed = await rtgtools_svdecompose(vcf, prefix="test", tbi=tbi)
    result = await rtgtools_bndeval(
        decomposed.vcf, vcf, prefix="test", query_tbi=decomposed.tbi, truth_tbi=tbi, regions_bed=regions
    )
    await check_md5s(
        {
            "tp": (result.tp_vcf, "976e8c05ba772bc4a4754921bdde77c3"),
            "tp tbi": (result.tp_tbi, EMPTY_TBI),
            "fn": (result.fn_vcf, "43692c2b0bf1ec3a7c4ed1deb1b13c47"),
            "fn tbi": (result.fn_tbi, EMPTY_TBI),
            "fp": (result.fp_vcf, "976e8c05ba772bc4a4754921bdde77c3"),
            "fp tbi": (result.fp_tbi, EMPTY_TBI),
            "baseline": (result.baseline_vcf, "43692c2b0bf1ec3a7c4ed1deb1b13c47"),
            "baseline tbi": (result.baseline_tbi, EMPTY_TBI),
            "weighted roc": (result.weighted_roc, "e66569caefbd0a1a9b01cd5c898a68c6"),
            "summary": (result.summary, "f4c8df93c8bdab603036bbc27b4a28c3"),
        }
    )


@env.task
async def test_rtgtools_cnveval_regions() -> None:
    query, query_tbi, truth, truth_tbi, regions = await asyncio.gather(
        fixture(SIMULATED_SV + ".vcf.gz"),
        fixture(SIMULATED_SV + ".vcf.gz.tbi"),
        fixture(SIMULATED_SV + "2.vcf.gz"),
        fixture(SIMULATED_SV + "2.vcf.gz.tbi"),
        fixture(HS + "genome/genome.bed"),
    )
    results = await rtgtools_cnveval_cmd.override(**TWO_CPUS)(
        query_vcf=query,
        query_tbi=query_tbi,
        truth_vcf=truth,
        truth_tbi=truth_tbi,
        evaluation_regions_bed=regions,
        prefix="test",
        args="",
    )
    result = await cnveval_result(results, "test")
    await check_md5s(
        {
            "baseline bed": (result.baseline_bed, "d93ffeef4c9bd38370f526980904dac6"),
            "calls bed": (result.calls_bed, "d93ffeef4c9bd38370f526980904dac6"),
            "weighted roc": (result.weighted_roc, "b9f4c2716ef1c8fb3341d1541982cfe9"),
            "summary": (result.summary, "f4c8df93c8bdab603036bbc27b4a28c3"),
        }
    )


tests = [
    test_rtgtools_format_fasta,
    test_rtgtools_format_fastqs,
    test_rtgtools_format_bam_rg,
    test_rtgtools_vcfeval_all_inputs,
    test_rtgtools_vcfeval_no_indexes,
    test_rtgtools_svdecompose_vcf,
    test_rtgtools_bndeval_regions,
    test_rtgtools_cnveval_regions,
]
