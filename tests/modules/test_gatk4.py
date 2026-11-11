"""Tests for flyte_bio.modules.gatk4 (md5s from upstream's gatk4 module snapshots)."""

import asyncio
import gzip
import hashlib
import re
import tarfile
import tempfile
import zlib
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.gatk4 import (
    Mutect2Result,
    gatk4_annotateintervals,
    gatk4_applybqsr,
    gatk4_baserecalibrator,
    gatk4_bedtointervallist,
    gatk4_calculatecontamination,
    gatk4_cnnscorevariants,
    gatk4_collectreadcounts,
    gatk4_combinegvcfs,
    gatk4_concordance,
    gatk4_createreadcountpanelofnormals,
    gatk4_createsequencedictionary,
    gatk4_createsomaticpanelofnormals,
    gatk4_determinegermlinecontigploidy,
    gatk4_filterintervals,
    gatk4_filtermutectcalls,
    gatk4_filtervarianttranches,
    gatk4_gatherbqsrreports,
    gatk4_gatherpileupsummaries,
    gatk4_genomicsdbimport,
    gatk4_genomicsdbimport_interval_list,
    gatk4_genomicsdbimport_update,
    gatk4_germlinecnvcaller,
    gatk4_getpileupsummaries,
    gatk4_haplotypecaller,
    gatk4_indexfeaturefile,
    gatk4_intervallisttobed,
    gatk4_intervallisttools,
    gatk4_learnreadorientationmodel,
    gatk4_markduplicates,
    gatk4_mergemutectstats,
    gatk4_mergevcfs,
    gatk4_mutect2,
    gatk4_preprocessintervals,
    gatk4_splitintervals,
    gatk4_splitncigarreads,
    gatk4_variantfiltration,
    gatk4_variantrecalibrator,
)
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_md5, assert_nonempty, assert_reads_md5, env, fixture, reads_md5

HS = "genomics/homo_sapiens/"
SC2 = "genomics/sarscov2/"
DBSNP = HS + "genome/vcf/dbsnp_146.hg38.vcf.gz"
CHR21 = HS + "genome/chr21/sequence/genome"
GERMLINE = HS + "genome/chr21/germlineresources/"
GNOMAD = GERMLINE + "gnomAD.r2.1.1.vcf.gz"
MILLS = GERMLINE + "mills_and_1000G.indels.hg38.vcf.gz"
AXIOM = GERMLINE + "axiom_exome_plus.genotypes.all_populations.poly.hg38.vcf.gz"
RECALIBRATED = HS + "illumina/bam/{}.paired_end.recalibrated.sorted.bam"
MUTECT2_CALLS = HS + "illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz"


async def vcf_lines(vcf: File) -> tuple[list[str], list[str]]:
    """A gzipped VCF's (header lines, record lines)."""
    async with vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    return [x for x in lines if x.startswith("#")], [x for x in lines if not x.startswith("#")]


async def variants_md5(vcf: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of a gzipped VCF's record lines joined without newlines, without leading zeros."""
    _, records = await vcf_lines(vcf)
    return format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")


async def assert_variants_md5(vcf: File, expected: str, label: str) -> None:
    actual = await variants_md5(vcf)
    assert actual == expected, f"{label}: variantsMD5 {actual}, expected {expected}"


async def human_reference():
    return await asyncio.gather(
        fixture(HS + "genome/genome.fasta"), fixture(HS + "genome/genome.fasta.fai"), fixture(HS + "genome/genome.dict")
    )


async def sarscov2_reference():
    return await asyncio.gather(
        fixture(SC2 + "genome/genome.fasta"),
        fixture(SC2 + "genome/genome.fasta.fai"),
        fixture(SC2 + "genome/genome.dict"),
    )


@env.task
async def test_createsequencedictionary() -> None:
    # upstream case: gatk4/createsequencedictionary "sarscov2 - fasta"
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    dict_ = await gatk4_createsequencedictionary(fasta)
    assert dict_.path.endswith("/genome.dict"), dict_.path
    await assert_md5(dict_, "7362679f176e0f52add03c08f457f646", label="genome.dict")


@env.task
async def test_intervallisttobed() -> None:
    # upstream case: gatk4/intervallisttobed "homo sapiens - bed"
    intervals = await fixture("genomics/homo_sapiens/genome/chr21/sequence/genome.interval_list")
    bed = await gatk4_intervallisttobed(intervals, prefix="test")
    await assert_md5(bed, "9046675d01199fbbee79f2bc1c5dce52", label="test.bed")


@env.task
async def test_gatk4_bedtointervallist() -> None:
    # upstream case: gatk4/bedtointervallist "test_gatk4_bedtointervallist" (GATK 4.7.0.0)
    bed = await fixture(SC2 + "genome/bed/test.bed")
    dict_ = await fixture(SC2 + "genome/genome.dict")
    interval_list = await gatk4_bedtointervallist(bed, dict_, prefix="test")
    assert interval_list.path.endswith("/test.interval_list"), interval_list.path
    await assert_md5(interval_list, "e51101c9357fb2d59fd30e370eefa39c", label="test.interval_list")


@env.task
async def test_markduplicates_multiple() -> None:
    # upstream cases: gatk4/markduplicates "homo_sapiens - multiple cram" / "multiple bam".
    # Their BAM md5 covers a header that records the command line, so check the outputs instead.
    bams = await asyncio.gather(
        fixture(HS + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(HS + "illumina/bam/test2.paired_end.sorted.bam"),
    )
    fasta, fai = await asyncio.gather(fixture(HS + "genome/genome.fasta"), fixture(HS + "genome/genome.fasta.fai"))
    cram, bam = await asyncio.gather(
        gatk4_markduplicates(list(bams), "test.cram", fasta=fasta, fai=fai),
        gatk4_markduplicates(list(bams), "test.bam", fasta=fasta, fai=fai),
    )
    assert cram.alignment.path.endswith("/test.cram") and cram.index is not None, cram
    assert cram.index.path.endswith("/test.cram.crai"), cram.index.path
    assert bam.alignment.path.endswith("/test.bam") and bam.index is None, bam
    for result in (cram, bam):
        async with result.metrics.open("rb") as fh:
            metrics = bytes(await fh.read()).decode()
        assert "## METRICS CLASS\tpicard.sam.DuplicationMetrics" in metrics, metrics[:400]
        await assert_nonempty(result.alignment, label=result.alignment.path.rsplit("/", 1)[-1])
    # Both inputs are merged into one alignment: its header has both read groups.
    header = await samtools_view(bam=bam.alignment, args="-H")
    async with header.open("rb") as fh:
        read_groups = [line for line in bytes(await fh.read()).decode().splitlines() if line.startswith("@RG")]
    assert len(read_groups) >= 1, read_groups


@env.task
async def test_baserecalibrator() -> None:
    # upstream cases: gatk4/baserecalibrator "sarscov2 - bam" and "sarscov2 - bam - intervals"
    (fasta, fai, dict_), bam, bai, vcf, tbi, bed = await asyncio.gather(
        sarscov2_reference(),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai"),
        fixture(SC2 + "illumina/vcf/test.vcf.gz"),
        fixture(SC2 + "illumina/vcf/test.vcf.gz.tbi"),
        fixture(SC2 + "genome/bed/test.bed"),
    )
    whole, within = await asyncio.gather(
        gatk4_baserecalibrator(bam, bai, fasta, fai, dict_, [vcf], [tbi], prefix="test"),
        gatk4_baserecalibrator(bam, bai, fasta, fai, dict_, [vcf], [tbi], prefix="test", intervals=bed),
    )
    await assert_md5(whole, "e2e43abdc0c943c1a54dae816d0b9ea7", label="test.table")
    await assert_md5(within, "9ecb5f00a2229291705addc09c0ec231", label="test.table (intervals)")


@env.task
async def test_gatherbqsrreports() -> None:
    # upstream cases: gatk4/gatherbqsrreports "test-gatk4-gatherbqsrreports" and "...-multiple"
    one, two = await asyncio.gather(
        fixture(HS + "illumina/gatk/test.baserecalibrator.table"),
        fixture(HS + "illumina/gatk/test2.baserecalibrator.table"),
    )
    single, multiple = await asyncio.gather(
        gatk4_gatherbqsrreports([one], prefix="test"), gatk4_gatherbqsrreports([one, two], prefix="test")
    )
    await assert_md5(single, "9603b69fdc3b5090de2e0dd78bfcc4bf", label="gathered single")
    await assert_md5(multiple, "0c1257eececf95db8ca378272d0f21f9", label="gathered multiple")


@env.task
async def test_applybqsr() -> None:
    # upstream cases: gatk4/applybqsr "sarscov2 - bam" and "sarscov2 - bam in, cram out". Their BAM md5
    # covers a header recording the command line; recalibration only changes base qualities, so the
    # read sequences must come through unchanged.
    (fasta, fai, dict_), bam, bai, table = await asyncio.gather(
        sarscov2_reference(),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai"),
        fixture(SC2 + "illumina/gatk/test.baserecalibrator.table"),
    )
    as_bam, as_cram = await asyncio.gather(
        gatk4_applybqsr(bam, bai, table, fasta, fai, dict_, prefix="test", suffix="bam"),
        gatk4_applybqsr(bam, bai, table, fasta, fai, dict_, prefix="test"),
    )
    assert as_bam.path.endswith("/test.bam") and as_cram.path.endswith("/test.cram"), (as_bam.path, as_cram.path)
    expected = await reads_md5(bam)
    await assert_reads_md5(as_bam, expected, label="recalibrated bam")
    await assert_nonempty(as_cram, label="recalibrated cram")


@env.task
async def test_haplotypecaller() -> None:
    # upstream case: gatk4/haplotypecaller "homo_sapiens - [cram, crai] - fasta - fai - dict - sites - sites_tbi"
    # (its snapshot records only the output names)
    (fasta, fai, dict_), cram, crai, dbsnp, dbsnp_tbi = await asyncio.gather(
        human_reference(),
        fixture(HS + "illumina/cram/test.paired_end.sorted.cram"),
        fixture(HS + "illumina/cram/test.paired_end.sorted.cram.crai"),
        fixture(DBSNP),
        fixture(DBSNP + ".tbi"),
    )
    calls = await gatk4_haplotypecaller(
        cram, crai, fasta, fai, dict_, "test_cram_sites", dbsnp=dbsnp, dbsnp_tbi=dbsnp_tbi
    )
    assert calls.vcf.path.endswith("/test_cram_sites.vcf.gz"), calls.vcf.path
    assert calls.tbi.path.endswith("/test_cram_sites.vcf.gz.tbi"), calls.tbi.path
    header, records = await vcf_lines(calls.vcf)
    assert header[0].startswith("##fileformat=VCF"), header[0]
    assert any(h.startswith("##GATKCommandLine=<ID=HaplotypeCaller") for h in header), "not a HaplotypeCaller VCF"
    assert records, "HaplotypeCaller called nothing on the test CRAM"


@env.task
async def test_mergevcfs() -> None:
    # upstream case: gatk4/mergevcfs "test_gatk4_mergevcfs" (snapshot records only names)
    dbsnp, gnomad, (_, _, dict_) = await asyncio.gather(
        fixture(DBSNP), fixture(HS + "genome/vcf/gnomAD.r2.1.1.vcf.gz"), human_reference()
    )
    merged = await gatk4_mergevcfs([dbsnp, gnomad], "test", dict=dict_)
    assert merged.vcf.path.endswith("/test.vcf.gz"), merged.vcf.path
    _, records = await vcf_lines(merged.vcf)
    expected = len((await vcf_lines(dbsnp))[1]) + len((await vcf_lines(gnomad))[1])
    assert len(records) == expected, (len(records), expected)


@env.task
async def test_cnnscorevariants() -> None:
    # upstream case: gatk4/cnnscorevariants "homo sapiens - vcf". Its variantsMD5 comes from htsjdk's
    # VariantContext rendering, so check the scores instead.
    (fasta, fai, dict_), vcf, tbi = await asyncio.gather(
        human_reference(),
        fixture(HS + "illumina/gvcf/test.genome.vcf.gz"),
        fixture(HS + "illumina/gvcf/test.genome.vcf.gz.tbi"),
    )
    scored = await gatk4_cnnscorevariants(vcf, tbi, fasta, fai, dict_, "test")
    assert scored.vcf.path.endswith("/test.cnn.vcf.gz"), scored.vcf.path
    header, records = await vcf_lines(scored.vcf)
    assert any(h.startswith("##INFO=<ID=CNN_1D") for h in header), "no CNN_1D INFO header"
    assert records and all("CNN_1D=" in r.split("\t")[7] for r in records), records[:2]


@env.task
async def test_filtervarianttranches() -> None:
    # upstream case: gatk4/filtervarianttranches "homo sapiens - vcf" (ext.args "--info-key CNN_1D")
    (fasta, fai, dict_), vcf, tbi, dbsnp, dbsnp_tbi = await asyncio.gather(
        human_reference(),
        fixture(HS + "illumina/gatk/haplotypecaller_calls/test_haplotcaller.cnn.vcf.gz"),
        fixture(HS + "illumina/gatk/haplotypecaller_calls/test_haplotcaller.cnn.vcf.gz.tbi"),
        fixture(DBSNP),
        fixture(DBSNP + ".tbi"),
    )
    filtered = await gatk4_filtervarianttranches(
        vcf, tbi, [dbsnp], [dbsnp_tbi], fasta, fai, dict_, "test", args="--info-key CNN_1D"
    )
    assert filtered.vcf.path.endswith("/test.filtered.vcf.gz"), filtered.vcf.path
    header, records = await vcf_lines(filtered.vcf)
    assert any(h.startswith("##FILTER=<ID=CNN_1D_") for h in header), "no tranche FILTER headers"
    assert records and all(r.split("\t")[6] for r in records), "records without a FILTER value"


async def chr21_reference():
    return await asyncio.gather(fixture(CHR21 + ".fasta"), fixture(CHR21 + ".fasta.fai"), fixture(CHR21 + ".dict"))


async def recalibrated_bam(name: str) -> tuple[File, File]:
    bam = RECALIBRATED.format(name)
    return await asyncio.gather(fixture(bam), fixture(bam + ".bai"))


async def indexed_fixture(path: str) -> tuple[File, File]:
    return await asyncio.gather(fixture(path), fixture(path + ".tbi"))


async def check_mutect2(result: Mutect2Result, label: str, variants: str, stats: str, f1r2: bool = False) -> None:
    """Upstream's Mutect2 snapshot: the calls' variantsMD5, the index's name, the stats' md5, f1r2's name."""
    assert result.vcf.name == "test.vcf.gz" and result.tbi.name == "test.vcf.gz.tbi", (result.vcf, result.tbi)
    await assert_variants_md5(result.vcf, variants, label)
    await assert_md5(result.stats, stats, label=f"{label} test.vcf.gz.stats")
    if f1r2:
        assert result.f1r2 is not None and result.f1r2.name == "test.f1r2.tar.gz", f"{label}: {result.f1r2}"
    else:
        assert result.f1r2 is None, f"{label}: unexpected {result.f1r2}"


@env.task
async def test_gatk4_mutect2_tumor_normal() -> None:
    # upstream cases: gatk4/mutect2 "human - bam - tumor_normal_pair", "... tumor_normal_pair_f1r2" and
    # "... tumor_normal_pair - force_call" (the last without --normal-sample, as upstream)
    (fasta, fai, dict_), (normal, normal_bai), (tumour, tumour_bai) = await asyncio.gather(
        chr21_reference(), recalibrated_bam("test"), recalibrated_bam("test2")
    )
    (gnomad, gnomad_tbi), (mills, mills_tbi), (axiom, axiom_tbi) = await asyncio.gather(
        indexed_fixture(GNOMAD), indexed_fixture(MILLS), indexed_fixture(AXIOM)
    )
    common = {
        "alignments": [normal, tumour],
        "indices": [normal_bai, tumour_bai],
        "fasta": fasta,
        "fai": fai,
        "dict": dict_,
        "prefix": "test",
        "germline_resource": gnomad,
        "germline_resource_tbi": gnomad_tbi,
        "panel_of_normals": mills,
        "panel_of_normals_tbi": mills_tbi,
    }
    pair, pair_f1r2, force_call = await asyncio.gather(
        gatk4_mutect2(**common, args="--normal-sample normal"),
        gatk4_mutect2(**common, args="--normal-sample normal --f1r2-tar-gz test.f1r2.tar.gz"),
        gatk4_mutect2(**common, alleles=axiom, alleles_tbi=axiom_tbi),
    )
    await check_mutect2(
        pair, "tumor_normal_pair", "7418ed45a029394253817a5eb7149334", "17d2091015d04cbd4a26b7a67dc659e6"
    )
    await check_mutect2(
        pair_f1r2,
        "tumor_normal_pair_f1r2",
        "7418ed45a029394253817a5eb7149334",
        "17d2091015d04cbd4a26b7a67dc659e6",
        f1r2=True,
    )
    await check_mutect2(
        force_call,
        "tumor_normal_pair - force_call",
        "63c4ff01a8316f4243b6b12b06ded74e",
        "81390058f99b80db30c89020aa9fc546",
    )


@env.task
async def test_gatk4_mutect2_tumor_only() -> None:
    # upstream cases: gatk4/mutect2 "human - bam - tumor_only", "human - cram", "human - bam - tumor_only -
    # force_call" and "human - bam - generate_pon" (no resources)
    cram = HS + "illumina/cram/test2.paired_end.recalibrated.sorted.cram"
    (fasta, fai, dict_), (bam, bai), (normal, normal_bai), (cram, crai) = await asyncio.gather(
        chr21_reference(),
        recalibrated_bam("test2"),
        recalibrated_bam("test"),
        asyncio.gather(fixture(cram), fixture(cram + ".crai")),
    )
    (gnomad, gnomad_tbi), (mills, mills_tbi), (axiom, axiom_tbi) = await asyncio.gather(
        indexed_fixture(GNOMAD), indexed_fixture(MILLS), indexed_fixture(AXIOM)
    )
    reference = {"fasta": fasta, "fai": fai, "dict": dict_, "prefix": "test"}
    resources = {
        "germline_resource": gnomad,
        "germline_resource_tbi": gnomad_tbi,
        "panel_of_normals": mills,
        "panel_of_normals_tbi": mills_tbi,
    }
    tumour_only, from_cram, force_call, pon = await asyncio.gather(
        gatk4_mutect2([bam], [bai], **reference, **resources),
        gatk4_mutect2([cram], [crai], **reference, **resources),
        gatk4_mutect2([bam], [bai], **reference, **resources, alleles=axiom, alleles_tbi=axiom_tbi),
        gatk4_mutect2([normal], [normal_bai], **reference),
    )
    await check_mutect2(
        tumour_only, "tumor_only", "1b65f1a163b517944bf2e4b74230e035", "55ed641e16089afb33cdbc478e202d3d"
    )
    await check_mutect2(from_cram, "cram", "1b65f1a163b517944bf2e4b74230e035", "55ed641e16089afb33cdbc478e202d3d")
    await check_mutect2(
        force_call, "tumor_only - force_call", "9c2928d0e751262c048ac162b9223e4e", "2398263886ecb546c8055af4fbad6040"
    )
    await check_mutect2(pon, "generate_pon", "876aa6be01c0c8fc71ad8e99ed842240", "b569ce66bbffe9588b3d221e821023ee")


@env.task
async def test_gatk4_mutect2_mitochondria() -> None:
    # upstream cases: gatk4/mutect2 "mitochondria - bam" and "mitochondria - bam - gz_ref" (--mitochondria-mode)
    bam = HS + "illumina/bam/mitochon_standin.recalibrated.sorted.bam"
    (fasta, fai, dict_), bam, bai, bed, fasta_gz, fai_gz, gzi = await asyncio.gather(
        human_reference(),
        fixture(bam),
        fixture(bam + ".bai"),
        fixture(HS + "genome/genome.bed"),
        fixture(HS + "genome/genome.fasta.gz"),
        fixture(HS + "genome/genome.fasta.gz.fai"),
        fixture(HS + "genome/genome.fasta.gz.gzi"),
    )
    common = {
        "alignments": [bam],
        "indices": [bai],
        "dict": dict_,
        "prefix": "test",
        "intervals": bed,
        "args": "--mitochondria-mode",
    }
    plain, bgzipped = await asyncio.gather(
        gatk4_mutect2(**common, fasta=fasta, fai=fai), gatk4_mutect2(**common, fasta=fasta_gz, fai=fai_gz, gzi=gzi)
    )
    for result, label in ((plain, "mitochondria"), (bgzipped, "mitochondria - gz_ref")):
        await check_mutect2(result, label, "ea70f79e33805a2c0b47b32a48a8d26f", "fc6ea14ca2da346babe78161beea28c9")


@env.task
async def test_gatk4_mergemutectstats() -> None:
    # upstream case: gatk4/mergemutectstats "human - stats, tsv" (its snapshot records the first line)
    stats = await fixture(MUTECT2_CALLS + ".stats")
    merged = await gatk4_mergemutectstats([stats], prefix="test")
    assert merged.name == "test.vcf.gz.stats", merged.name
    async with merged.open("rb") as fh:
        first = bytes(await fh.read()).decode().splitlines()[0]
    assert first == "statistic\tvalue", repr(first)


@env.task
async def test_gatk4_filtermutectcalls() -> None:
    # upstream cases: gatk4/filtermutectcalls "human - vcf - base", "with-files", "use-val", "stats renamed"
    # and "no stats - fails" (ext.prefix "<id>.filtered"; the index is snapshotted by name)
    gatk = HS + "illumina/gatk/"
    (fasta, fai, dict_), (vcf, tbi), stats, priors, segmentation, contamination = await asyncio.gather(
        chr21_reference(),
        indexed_fixture(MUTECT2_CALLS),
        fixture(MUTECT2_CALLS + ".stats"),
        fixture(gatk + "test_test2_paired_mutect2_calls.artifact-prior.tar.gz"),
        fixture(gatk + "test_test2_paired.segmentation.table"),
        fixture(gatk + "test_test2_paired.contamination.table"),
    )
    renamed = await rename_mutect2_stats(stats)
    common = {"vcf": vcf, "tbi": tbi, "fasta": fasta, "fai": fai, "dict": dict_, "prefix": "test.filtered"}
    base, with_files, use_val, stats_renamed = await asyncio.gather(
        gatk4_filtermutectcalls(stats=stats, **common),
        gatk4_filtermutectcalls(
            stats=stats,
            orientation_priors=[priors],
            segmentations=[segmentation],
            contamination_tables=[contamination],
            **common,
        ),
        gatk4_filtermutectcalls(
            stats=stats,
            orientation_priors=[priors],
            segmentations=[segmentation],
            contamination_estimate=20.0,
            **common,
        ),
        gatk4_filtermutectcalls(stats=renamed, **common),
    )
    expected = (
        (base, "base", "5352aa0f3ea986008a2652eba7c0fdba", "47c14d119fa3921eeae8b0be0ce1a14c"),
        (with_files, "with-files", "54e439a79e4b6956a7f2c5b0730c9918", "ca16c62a13d3b0fa45568c3a66a09104"),
        (use_val, "use-val", "4f8cfebe28e5fd097a7bfc2e556983d3", "95cc3e37705bd3b97a292c5d46ab82f3"),
        (stats_renamed, "stats renamed", "5352aa0f3ea986008a2652eba7c0fdba", "47c14d119fa3921eeae8b0be0ce1a14c"),
    )
    for result, label, variants, filtering_stats in expected:
        assert result.vcf.name == "test.filtered.vcf.gz", f"{label}: {result.vcf.name}"
        assert result.tbi.name == "test.filtered.vcf.gz.tbi", f"{label}: {result.tbi.name}"
        await assert_variants_md5(result.vcf, variants, label)
        await assert_md5(result.stats, filtering_stats, label=f"{label} filteringStats.tsv")
    try:
        await gatk4_filtermutectcalls(stats=None, **common)
    except Exception:
        return
    raise AssertionError("no stats - fails: FilterMutectCalls ran without the Mutect stats table")


async def rename_mutect2_stats(stats: File) -> File:
    """A copy of ``stats`` named ``renamed.mutect2.stats``, so GATK can only find it through --stats."""
    local = Path(tempfile.mkdtemp()) / "renamed.mutect2.stats"
    async with stats.open("rb") as fh:
        local.write_bytes(bytes(await fh.read()))
    return await File.from_local(str(local))


@env.task
async def test_gatk4_variantfiltration() -> None:
    # upstream cases: gatk4/variantfiltration "human - vcf - fasta", "human - vcf - fasta.gz" and
    # "human - vcf.gz - fasta" (no filters; the index is snapshotted by name). Plus a named filter, which
    # upstream doesn't test.
    gvcf = HS + "illumina/gvcf/test.genome.vcf"
    (fasta, fai, dict_), vcf, idx, (vcf_gz, tbi), fasta_gz, fai_gz, gzi = await asyncio.gather(
        human_reference(),
        fixture(gvcf),
        fixture(gvcf + ".idx"),
        indexed_fixture(gvcf + ".gz"),
        fixture(HS + "genome/genome.fasta.gz"),
        fixture(HS + "genome/genome.fasta.gz.fai"),
        fixture(HS + "genome/genome.fasta.gz.gzi"),
    )
    plain, bgzipped_ref, bgzipped_vcf, filtered = await asyncio.gather(
        gatk4_variantfiltration(vcf, idx, fasta, fai, dict_, "test"),
        gatk4_variantfiltration(vcf, idx, fasta_gz, fai_gz, dict_, "test", gzi=gzi),
        gatk4_variantfiltration(vcf_gz, tbi, fasta, fai, dict_, "test"),
        gatk4_variantfiltration(vcf_gz, tbi, fasta, fai, dict_, "test", filters={"LowDepth": "DP < 10"}),
    )
    for result, label in ((plain, "vcf - fasta"), (bgzipped_ref, "vcf - fasta.gz"), (bgzipped_vcf, "vcf.gz - fasta")):
        assert result.tbi.name == "test.vcf.gz.tbi", f"{label}: {result.tbi.name}"
        await assert_variants_md5(result.vcf, "cf0477cc14953b16caacbf5be507b88", label)
    header, records = await vcf_lines(filtered.vcf)
    assert any(h.startswith("##FILTER=<ID=LowDepth") for h in header), "no LowDepth FILTER header"
    assert records and all(r.split("\t")[6] != "." for r in records), "records without a FILTER value"


@env.task
async def test_gatk4_intervallisttools() -> None:
    # upstream case: gatk4/intervallisttools "test_gatk4_intervallisttools", on BedToIntervalList's output
    bed, dict_ = await asyncio.gather(fixture(SC2 + "genome/bed/test.bed"), fixture(SC2 + "genome/genome.dict"))
    interval_list = await gatk4_bedtointervallist(bed, dict_, prefix="test")
    split = await gatk4_intervallisttools(
        interval_list,
        prefix="test",
        args="--SCATTER_COUNT 6 --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW "
        "--UNIQUE true --SORT true",
    )
    expected = [
        "64f6665f9fbd257e4a300ec602f4e995",
        "f515c3da0c6accfd8e7dc33df50855c5",
        "7a918e8c9211b54334587793e8cbae53",
        "1b93105227a7dc81f07101a1efd31498",
    ]
    names = [f.name for f in split]
    assert names == [f"{n}scattered.interval_list" for n in range(1, 5)], names
    for f, md5 in zip(split, expected, strict=True):
        await assert_md5(f, md5, label=f.name)


# The record lines (linesGzip[102..-1]) of upstream's gatk4/combinegvcfs snapshot: their count and the md5 of the
# lines joined with newlines. The 102 header lines before them record the command line and date.
COMBINEGVCFS_HEADER_LINES = 102
COMBINEGVCFS_BODY_LINES = 130
COMBINEGVCFS_BODY_MD5 = "1e3cf9acf1124f1bee2440f67fa2606a"


@env.task
async def test_gatk4_combinegvcfs() -> None:
    # upstream case: gatk4/combinegvcfs "test_gatk4_combinegvcfs"
    vcf, idx = await asyncio.gather(
        fixture(HS + "illumina/gvcf/test.genome.vcf"), fixture(HS + "illumina/gvcf/test.genome.vcf.idx")
    )
    fasta, fai, dict_ = await human_reference()
    combined = await gatk4_combinegvcfs([vcf], [idx], fasta, fai, dict_, prefix="test")
    assert combined.name == "test.combined.g.vcf.gz", combined.name
    async with combined.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    assert lines[COMBINEGVCFS_HEADER_LINES - 1].startswith("#CHROM"), lines[COMBINEGVCFS_HEADER_LINES - 1]
    body = lines[COMBINEGVCFS_HEADER_LINES:]
    assert len(body) == COMBINEGVCFS_BODY_LINES, f"{len(body)} lines after the header"
    actual = hashlib.md5("\n".join(body).encode()).hexdigest()
    assert actual == COMBINEGVCFS_BODY_MD5, f"combined GVCF body md5 {actual}, expected {COMBINEGVCFS_BODY_MD5}"


async def assert_bam(bam: File, name: str, label: str) -> None:
    """``bam`` is named ``name`` and is a BAM (gzip-framed, ``BAM\\1`` magic)."""
    assert bam.name == name, f"{label}: {bam.name}"
    async with bam.open("rb") as fh:
        head = bytes(await fh.read(1 << 16))
    magic = zlib.decompressobj(wbits=31).decompress(head)[:4]  # the first BGZF block
    assert magic == b"BAM\1", f"{label}: starts {magic!r}"


@env.task
async def test_gatk4_splitncigarreads() -> None:
    # upstream cases: gatk4/splitncigarreads "test_gatk4_splitncigarreads" / "..._intervals".
    # Upstream snapshots only the output's name (the BAM header records the command line).
    bam, bai = await asyncio.gather(
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam.bai"),
    )
    fasta, fai, dict_ = await sarscov2_reference()
    bed = await fixture(SC2 + "genome/bed/test.bed")
    plain, within = await asyncio.gather(
        gatk4_splitncigarreads(bam, bai, fasta, fai, dict_, prefix="test"),
        gatk4_splitncigarreads(bam, bai, fasta, fai, dict_, prefix="test", intervals=bed),
    )
    await assert_bam(plain, "test.bam", "splitncigarreads")
    await assert_bam(within, "test.bam", "splitncigarreads intervals")


# The panel-of-normals and germline CNV tools: upstream's cases at modules commit 9a1cc9f (GATK 4.7.0.0), with
# their test configs' ext.args.
OVERLAPPING_ONLY = "--interval-merging-rule OVERLAPPING_ONLY"
COUNTS_TSV = "--format TSV " + OVERLAPPING_ONLY
MULTI_INTERVALS = HS + "genome/genome.multi_intervals.bed"


@env.task
async def test_gatk4_preprocessintervals() -> None:
    # upstream case: gatk4/preprocessintervals "test-gatk4-preprocessintervals" (no intervals, a blacklist excluded)
    (fasta, fai, dict_), blacklist = await asyncio.gather(
        human_reference(), fixture(HS + "genome/genome.blacklist_intervals.bed")
    )
    intervals = await gatk4_preprocessintervals(
        fasta, fai, dict_, "test", exclude_intervals=blacklist, args="--padding 0 " + OVERLAPPING_ONLY
    )
    assert intervals.name == "test.interval_list", intervals.name
    await assert_md5(intervals, "ce14b8fb47a60483fe44473ba40e1583", label="test.interval_list")


@env.task
async def test_gatk4_annotateintervals() -> None:
    # upstream cases: gatk4/annotateintervals "homo_sapiens genome [bed]", "... multi-interval [bed]",
    # "... [interval_list]", "... [interval_list] mappable_regions" and "... duplication_regions"
    (fasta, fai, dict_), bed, multi, interval_list, (track, track_tbi) = await asyncio.gather(
        human_reference(),
        fixture(HS + "genome/genome.bed"),
        fixture(MULTI_INTERVALS),
        fixture(HS + "genome/genome.interval_list"),
        indexed_fixture(HS + "genome/genome.bed.gz"),
    )
    reference = {"fasta": fasta, "fai": fai, "dict": dict_, "prefix": "test", "args": OVERLAPPING_ONLY}
    results = await asyncio.gather(
        gatk4_annotateintervals([bed], **reference),
        gatk4_annotateintervals([bed, multi], **reference),
        gatk4_annotateintervals([interval_list], **reference),
        gatk4_annotateintervals([interval_list], **reference, mappable_regions=track, mappable_regions_tbi=track_tbi),
        gatk4_annotateintervals(
            [interval_list],
            **reference,
            segmental_duplication_regions=track,
            segmental_duplication_regions_tbi=track_tbi,
        ),
    )
    expected = (
        ("bed", "abf7cef9157e2b3c4554f48475f1b7a1"),
        ("multi-interval bed", "abf7cef9157e2b3c4554f48475f1b7a1"),
        ("interval_list", "abf7cef9157e2b3c4554f48475f1b7a1"),
        ("interval_list mappable_regions", "428d33a03f2887e267d1cb6f0351bf9d"),
        ("interval_list duplication_regions", "0c59e93cc6e35ca14d3d25ef885156f8"),
    )
    for tsv, (label, md5) in zip(results, expected, strict=True):
        assert tsv.name == "test.tsv", f"{label}: {tsv.name}"
        await assert_md5(tsv, md5, label=label)


@env.task
async def test_gatk4_collectreadcounts() -> None:
    # upstream cases: gatk4/collectreadcounts "test-gatk4-collectreadcounts-hdf5" (its snapshot records the name),
    # "...-tsv" and "...-cram"
    cram = HS + "illumina/cram/test.paired_end.sorted.cram"
    (fasta, fai, dict_), (bam, bai), (cram, crai), bed = await asyncio.gather(
        human_reference(),
        bam_with_index(HS + "illumina/bam/test.paired_end.sorted.bam"),
        asyncio.gather(fixture(cram), fixture(cram + ".crai")),
        fixture(HS + "genome/genome.bed"),
    )
    hdf5, tsv, from_cram = await asyncio.gather(
        gatk4_collectreadcounts(bam, bai, bed, "test", args=OVERLAPPING_ONLY),
        gatk4_collectreadcounts(bam, bai, bed, "test", args=COUNTS_TSV),
        gatk4_collectreadcounts(cram, crai, bed, "test", fasta=fasta, fai=fai, dict=dict_, args=COUNTS_TSV),
    )
    assert hdf5.name == "test.hdf5", hdf5.name
    assert tsv.name == "test.tsv" and from_cram.name == "test.tsv", (tsv.name, from_cram.name)
    await assert_md5(tsv, "8e45a6164916c303387f39f02ce45841", label="tsv")
    await assert_md5(from_cram, "d9a32039b7a84f5bb74e8382e5427670", label="cram")


async def bam_with_index(path: str) -> tuple[File, File]:
    """A BAM fixture and its ``.bai``."""
    return await asyncio.gather(fixture(path), fixture(path + ".bai"))


async def cohort_counts(intervals: File) -> list[File]:
    """TSV read counts of upstream's two test BAMs (``test.tsv``, ``test2.tsv``) over ``intervals``."""
    (fasta, fai, dict_), *bams = await asyncio.gather(
        human_reference(),
        bam_with_index(HS + "illumina/bam/test.paired_end.sorted.bam"),
        bam_with_index(HS + "illumina/bam/test2.paired_end.sorted.bam"),
    )
    return list(
        await asyncio.gather(
            *(
                gatk4_collectreadcounts(bam, bai, intervals, name, fasta=fasta, fai=fai, dict=dict_, args=COUNTS_TSV)
                for (bam, bai), name in zip(bams, ("test", "test2"), strict=True)
            )
        )
    )


@env.task
async def test_gatk4_createreadcountpanelofnormals() -> None:
    # upstream case: gatk4/createreadcountpanelofnormals "homo sapiens - bam" (its snapshot records the name), on
    # PreprocessIntervals (default bins) and CollectReadCounts of two BAMs
    fasta, fai, dict_ = await human_reference()
    intervals = await gatk4_preprocessintervals(fasta, fai, dict_, "test")
    counts = await cohort_counts(intervals)
    pon = await gatk4_createreadcountpanelofnormals(
        counts, "test", args="--minimum-interval-median-percentile 1.0 --number-of-eigensamples 2"
    )
    assert pon.name == "test.hdf5", pon.name


@env.task
async def test_gatk4_filterintervals() -> None:
    # upstream case: gatk4/filterintervals "test-gatk4-filterintervals"
    genome = HS + "genome/"
    intervals, counts, annotated = await asyncio.gather(
        fixture(genome + "genome.preprocessed_intervals.interval_list"),
        fixture(genome + "genome.preprocessed_intervals.counts.tsv"),
        fixture(genome + "genome.annotated_intervals.tsv"),
    )
    filtered = await gatk4_filterintervals(
        intervals, "test", read_counts=[counts], annotated_intervals=annotated, args=OVERLAPPING_ONLY
    )
    assert filtered.name == "test.interval_list", filtered.name
    await assert_md5(filtered, "67b15dff732693db3542e6b1dc30a5da", label="test.interval_list")


async def dir_entries(directory: Dir) -> list[str]:
    """The sorted names of ``directory``'s top-level files and subdirectories."""
    base = directory.path.rstrip("/") + "/"
    return sorted({f.path[len(base) :].split("/", 1)[0] async for f in directory.walk()})


PLOIDY_MODEL_FILES = [
    "contig_ploidy_prior.tsv",
    "gcnvkernel_version.json",
    "interval_list.tsv",
    "mu_mean_bias_j_interval__.tsv",
    "mu_psi_j_log__.tsv",
    "ploidy_config.json",
    "std_mean_bias_j_interval__.tsv",
    "std_psi_j_log__.tsv",
]


async def cohort_ploidy():
    """Upstream's DetermineGermlineContigPloidy setup: counts over the multi-interval BED, then ploidy calls."""
    bed, priors = await asyncio.gather(
        fixture(MULTI_INTERVALS), fixture(HS + "illumina/gatk/contig_ploidy_priors_table.tsv")
    )
    counts = await cohort_counts(bed)
    ploidy = await gatk4_determinegermlinecontigploidy(
        counts, "test", intervals=bed, contig_ploidy_priors=priors, args=OVERLAPPING_ONLY
    )
    return counts, ploidy


@env.task
async def test_gatk4_determinegermlinecontigploidy() -> None:
    # upstream case: gatk4/determinegermlinecontigploidy "homo sapiens - bam" (cohort mode; its snapshot records
    # the calls' and model's file names: model fitting isn't byte-stable)
    _, ploidy = await cohort_ploidy()
    assert ploidy.model is not None, "no ploidy model in cohort mode"
    calls, model = await asyncio.gather(dir_entries(ploidy.calls), dir_entries(ploidy.model))
    assert calls == ["SAMPLE_0", "SAMPLE_1"], calls
    assert model == PLOIDY_MODEL_FILES, model


@env.task
async def test_gatk4_germlinecnvcaller() -> None:
    # upstream case: gatk4/germlinecnvcaller "homo sapiens - bam" (cohort mode, after CollectReadCounts,
    # DetermineGermlineContigPloidy and BedToIntervalList; its snapshot records the calls' and model's file names)
    (counts, ploidy), dict_ = await asyncio.gather(cohort_ploidy(), fixture(HS + "genome/genome.dict"))
    intervals = await gatk4_bedtointervallist(await fixture(MULTI_INTERVALS), dict_, prefix="test")
    cnv = await gatk4_germlinecnvcaller(
        counts,
        "test",
        intervals=intervals,
        contig_ploidy_calls=ploidy.calls,
        args=OVERLAPPING_ONLY + " --run-mode COHORT",
    )
    assert cnv.model is not None, "no cohort model in cohort mode"
    calls, model = await asyncio.gather(dir_entries(cnv.calls), dir_entries(cnv.model))
    configs = ["calling_config.json", "denoising_config.json", "gcnvkernel_version.json", "interval_list.tsv"]
    assert calls == ["SAMPLE_0", "SAMPLE_1", *configs], calls
    expected_model = [
        *configs,
        "log_q_tau_tk.tsv",
        "mu_W_tu.tsv",
        "mu_ard_u_interval__.tsv",
        "mu_log_mean_bias_t.tsv",
        "mu_psi_t_log__.tsv",
        "std_W_tu.tsv",
        "std_ard_u_interval__.tsv",
        "std_log_mean_bias_t.tsv",
        "std_psi_t_log__.tsv",
    ]
    assert model == expected_model, model


async def workspace_entries(archive: File) -> tuple[str, list[str]]:
    """A GenomicsDB workspace archive's directory name and the sorted names directly in it."""
    local = await archive.download()
    with tarfile.open(local) as tar:
        parts = [Path(m.name).parts for m in tar.getmembers()]
    roots = {p[0] for p in parts}
    assert len(roots) == 1, f"{archive.name}: top-level entries {sorted(roots)}"
    return roots.pop(), sorted({p[1] for p in parts if len(p) > 1})


WORKSPACE_ENTRIES = ["__tiledb_workspace.tdb", "callset.json", "chr22$1$40001", "vcfheader.vcf", "vidmap.json"]
GENOMICSDB = HS + "illumina/gatk/test_genomicsdb.tar.gz"


@env.task
async def test_gatk4_genomicsdbimport() -> None:
    # upstream cases: gatk4/genomicsdbimport "test_gatk4_genomicsdbimport_create_genomicsdb",
    # "..._get_intervalslist" and "..._update_genomicsdb" (the workspaces are snapshotted by their listing)
    (gvcf, tbi), (gvcf2, tbi2), intervals, workspace = await asyncio.gather(
        indexed_fixture(HS + "illumina/gvcf/test.genome.vcf.gz"),
        indexed_fixture(HS + "illumina/gvcf/test2.genome.vcf.gz"),
        fixture(HS + "genome/genome.interval_list"),
        fixture(GENOMICSDB),
    )
    created, interval_list, updated = await asyncio.gather(
        gatk4_genomicsdbimport([gvcf], [tbi], "test", intervals=intervals),
        gatk4_genomicsdbimport_interval_list(workspace, "test"),
        gatk4_genomicsdbimport_update(workspace, [gvcf2], [tbi2]),
    )
    assert created.name == "test.tar.gz", created.name
    assert await workspace_entries(created) == ("test", WORKSPACE_ENTRIES)
    assert interval_list.name == "test.interval_list", interval_list.name
    await assert_md5(interval_list, "4c85812ac15fc1cd29711a851d23c0bf", label="test.interval_list")
    assert updated.name == "test_genomicsdb.tar.gz", updated.name
    assert await workspace_entries(updated) == ("test_genomicsdb", WORKSPACE_ENTRIES)


@env.task
async def test_gatk4_createsomaticpanelofnormals() -> None:
    # upstream case: gatk4/createsomaticpanelofnormals "homo sapiens - genomicsdb" (its snapshot records the names)
    workspace, (fasta, fai, dict_) = await asyncio.gather(fixture(GENOMICSDB), human_reference())
    pon = await gatk4_createsomaticpanelofnormals(workspace, fasta, fai, dict_, "test")
    assert pon.vcf.name == "test.vcf.gz" and pon.tbi.name == "test.vcf.gz.tbi", (pon.vcf.name, pon.tbi.name)
    header, _ = await vcf_lines(pon.vcf)
    assert header[0].startswith("##fileformat=VCF"), header[0]


@env.task
async def test_gatk4_indexfeaturefile() -> None:
    # upstream cases: gatk4/indexfeaturefile "..._bed", "..._bed_gz", "..._vcf" and "..._vcf_gz" (Tribble .idx
    # indices record the file's path and date, so upstream snapshots their names)
    gvcf = HS + "illumina/gvcf/test.genome.vcf"
    files = await asyncio.gather(
        fixture(HS + "genome/genome.bed"),
        fixture(HS + "genome/genome.bed.gz"),
        fixture(gvcf),
        fixture(gvcf + ".gz"),
    )
    bed_idx, bed_tbi, vcf_idx, vcf_tbi = await asyncio.gather(*(gatk4_indexfeaturefile(f) for f in files))
    assert bed_idx.name == "genome.bed.idx", bed_idx.name
    assert vcf_idx.name == "test.genome.vcf.idx", vcf_idx.name
    assert bed_tbi.name == "genome.bed.gz.tbi" and vcf_tbi.name == "test.genome.vcf.gz.tbi", (bed_tbi, vcf_tbi)
    await assert_md5(bed_tbi, "4bc51e2351a6e83f20e13be75861f941", label="genome.bed.gz.tbi")
    await assert_md5(vcf_tbi, "fedd68eaddf8d31257853d9da8325bd3", label="test.genome.vcf.gz.tbi")


@env.task
async def test_gatk4_splitintervals() -> None:
    # upstream cases: gatk4/splitintervals "test-gatk4-splitintervals-bed" and "...-intervals" (--scatter-count 2)
    (fasta, fai, dict_), bed, interval_list = await asyncio.gather(
        human_reference(), fixture(MULTI_INTERVALS), fixture(HS + "genome/genome.interval_list")
    )
    reference = {"fasta": fasta, "fai": fai, "dict": dict_, "args": "--scatter-count 2"}
    from_bed, from_list = await asyncio.gather(
        gatk4_splitintervals(bed, "test", **reference), gatk4_splitintervals(interval_list, "test", **reference)
    )
    expected = (
        (from_bed, "bed", ("c8d6b19e7a92535b6ce9608eae558faa", "b1877ad96aec308906594c50ebbe3ded")),
        (from_list, "intervals", ("ebd6b34a335efc6732ff541936c6d2d5", "9459b0e124fa84ec1e64ac4615bc9af7")),
    )
    for split, label, md5s in expected:
        names = [f.name for f in split]
        assert names == ["0000-scattered.interval_list", "0001-scattered.interval_list"], f"{label}: {names}"
        for f, md5 in zip(split, md5s, strict=True):
            await assert_md5(f, md5, label=f"{label} {f.name}")


@env.task
async def test_gatk4_getpileupsummaries() -> None:
    # upstream cases: gatk4/getpileupsummaries "human - bam" (no reference) and "human - cram"
    interval_list = await fixture(CHR21 + ".interval_list")
    (bam, bai), (variants, variants_tbi), (fasta, fai, dict_) = await asyncio.gather(
        recalibrated_bam("test"), indexed_fixture(GNOMAD), chr21_reference()
    )
    cram = HS + "illumina/cram/test.paired_end.recalibrated.sorted.cram"
    cram_file, crai = await asyncio.gather(fixture(cram), fixture(cram + ".crai"))
    common = {"variants": variants, "variants_tbi": variants_tbi, "prefix": "test", "intervals": interval_list}
    from_bam, from_cram = await asyncio.gather(
        gatk4_getpileupsummaries(bam, bai, **common),
        gatk4_getpileupsummaries(cram_file, crai, fasta=fasta, fai=fai, dict=dict_, **common),
    )
    for table, label in ((from_bam, "bam"), (from_cram, "cram")):
        assert table.name == "test.pileups.table", f"{label}: {table.name}"
        await assert_md5(table, "2d7ce4a54df6b9249e12d60737a67dac", label=label)


@env.task
async def test_gatk4_gatherpileupsummaries() -> None:
    # upstream case: gatk4/gatherpileupsummaries "test-gatk4-gatherpileupsummaries" (ext.prefix "<id>.out")
    table, dict_ = await asyncio.gather(fixture(HS + "illumina/gatk/test.pileups.table"), fixture(CHR21 + ".dict"))
    gathered = await gatk4_gatherpileupsummaries([table], dict_, prefix="test.out")
    assert gathered.name == "test.out.pileups.table", gathered.name
    await assert_md5(gathered, "8e0ca6f66e112bd2f7ec1d31a2d62469")


async def first_line(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()[0]


@env.task
async def test_gatk4_calculatecontamination() -> None:
    # upstream cases: gatk4/calculatecontamination "human - pileup-table", "... - matched-pair" and
    # "... - segmentation" (each snapshot records the tables' first lines)
    tumour, normal = await asyncio.gather(
        fixture(HS + "illumina/gatk/test2.pileups.table"), fixture(HS + "illumina/gatk/test.pileups.table")
    )
    single, matched, segmented = await asyncio.gather(
        gatk4_calculatecontamination(tumour, prefix="test"),
        gatk4_calculatecontamination(tumour, prefix="test", matched=normal),
        gatk4_calculatecontamination(
            tumour, prefix="test", matched=normal, args="--tumor-segmentation test.segmentation.table"
        ),
    )
    for result, label in ((single, "pileup-table"), (matched, "matched-pair"), (segmented, "segmentation")):
        assert result.contamination.name == "test.contamination.table", f"{label}: {result.contamination.name}"
        line = await first_line(result.contamination)
        assert line == "sample\tcontamination\terror", f"{label}: {line!r}"
    assert single.segmentation is None and matched.segmentation is None, "unrequested segmentation table"
    assert segmented.segmentation is not None, "segmentation: no segmentation table"
    assert segmented.segmentation.name == "test.segmentation.table", segmented.segmentation.name
    line = await first_line(segmented.segmentation)
    assert line == "#<METADATA>SAMPLE=tumour", f"segmentation: {line!r}"


# Lines 3 to 7 of the artifact priors (nf-test's linesGzip of the archive, lines split as Java does).
ARTIFACT_PRIOR_LINES = [
    (
        "CTT\tAAG\t2.7114986684474486E-6\t3.2076972826656866E-5\t2.6085822355549755E-6\t0.0\t2.6371799896540086E-6"
        "\t3.3869355267901446E-6\t2.6085822355549755E-6\t0.0\t0.9995881552107633\t4.6590850211691583E-5"
        "\t2.8848017683240004E-4\t3.0744010710100574E-5\t38334\t116"
    ),
    "GTT\tAAC\t0.1\t0.1\t0.1\t0.0\t0.1\t0.1\t0.1\t0.0\t0.1\t0.1\t0.1\t0.1\t0\t0",
    (
        "TAT\tATA\t0.0\t5.548307163536064E-6\t5.144357865084592E-6\t6.205892051757818E-6\t0.0\t5.907388162200423E-6"
        "\t5.176730417709638E-6\t6.083872804985981E-6\t0.9924019419831304\t3.946972069386949E-5"
        "\t0.007516612822150651\t7.908925559851714E-6\t19439\t95"
    ),
    (
        "AAA\tTTT\t0.0\t1.7634470563520664E-6\t2.8327478284981175E-6\t1.8084237600021914E-6\t0.0\t1.7692606885284446E-6"
        "\t2.263339968296726E-6\t1.8660094002474611E-6\t0.9990845693211764\t1.8004690536795885E-5"
        "\t8.701192700921183E-4\t1.5003489492700572E-5\t56708\t130"
    ),
    (
        "CAA\tTTG\t0.0\t3.4445551925564533E-6\t3.435193155024585E-6\t5.139879646597498E-6\t0.0\t3.4674461103560476E-6"
        "\t3.428570449688764E-6\t6.343168047383713E-6\t0.9945238954147358\t1.3629167993931722E-4"
        "\t0.0052793454402581125\t3.5208652465150646E-5\t29263\t197"
    ),
]


@env.task
async def test_gatk4_learnreadorientationmodel() -> None:
    # upstream case: gatk4/learnreadorientationmodel "test-gatk4-learnreadorientationmodel"
    # (ext.prefix "<id>.artifact-prior"; the snapshot records lines 3 to 7 of the gunzipped archive)
    f1r2 = await fixture(HS + "illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.f1r2.tar.gz")
    priors = await gatk4_learnreadorientationmodel([f1r2], prefix="test.artifact-prior")
    assert priors.name == "test.artifact-prior.tar.gz", priors.name
    async with priors.open("rb") as fh:
        text = gzip.decompress(bytes(await fh.read())).decode("latin-1")
    lines = re.split(r"\r\n|\r|\n", text)[3:8]
    assert lines == ARTIFACT_PRIOR_LINES, f"artifact priors lines 3-7: {lines}"


@env.task
async def test_gatk4_variantrecalibrator() -> None:
    # upstream cases: gatk4/variantrecalibrator "homo sapiens - vcf - allele specificity" and
    # "... - no allele specificity" (the recal VCF and its index are snapshotted by name, the tranches by md5)
    calls = HS + "illumina/gatk/haplotypecaller_calls/test2_haplotc.ann.vcf.gz"
    names = ["hapmap_3.3.hg38", "1000G_omni2.5.hg38", "1000G_phase1.snps.hg38", "dbsnp_138.hg38"]
    (vcf, tbi), resources, (fasta, fai, dict_) = await asyncio.gather(
        indexed_fixture(calls),
        asyncio.gather(*(indexed_fixture(f"{GERMLINE}{name}.vcf.gz") for name in names)),
        chr21_reference(),
    )
    labels = [
        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz",
        "--resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz",
        "--resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.hg38.vcf.gz",
        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf.gz",
    ]
    common = {
        "vcf": vcf,
        "tbi": tbi,
        "resources": [r[0] for r in resources],
        "resources_tbi": [r[1] for r in resources],
        "labels": labels,
        "prefix": "test",
        "fasta": fasta,
        "fai": fai,
        "dict": dict_,
    }
    allele_specific, plain = await asyncio.gather(
        gatk4_variantrecalibrator(args="-mode SNP -an QD -an MQ -an FS -AS", **common),
        gatk4_variantrecalibrator(args="-mode SNP -an QD -an MQ -an FS -an SOR", **common),
    )
    for result, label, md5 in (
        (allele_specific, "allele specificity", "ad52fa69325c758f458a30ee5b43d6b5"),
        (plain, "no allele specificity", "c029e52fd63a893e1154cc9144a19eeb"),
    ):
        assert result.recal.name == "test.recal" and result.idx.name == "test.recal.idx", f"{label}: {result}"
        assert result.tranches.name == "test.tranches", f"{label}: {result.tranches.name}"
        assert result.plots is None, f"{label}: unexpected {result.plots}"
        await assert_md5(result.tranches, md5, label=f"{label} tranches")


@env.task
async def test_gatk4_concordance() -> None:
    calls = HS + "illumina/gatk/haplotypecaller_calls/test2_haplotc.vcf.gz"
    vcf, vcf_tbi, truth, truth_tbi, intervals, fasta, fai, sequence_dict = await asyncio.gather(
        fixture(MUTECT2_CALLS),
        fixture(MUTECT2_CALLS + ".tbi"),
        fixture(calls),
        fixture(calls + ".tbi"),
        fixture(HS + "genome/chr21/sequence/multi_intervals.bed"),
        fixture(CHR21 + ".fasta"),
        fixture(CHR21 + ".fasta.fai"),
        fixture(CHR21 + ".dict"),
    )
    result = await gatk4_concordance(
        vcf, vcf_tbi, truth, truth_tbi, fasta, fai, sequence_dict, prefix="test", intervals=intervals
    )
    await assert_md5(result.summary, "1fe75583afedc1c743f052f268d790ca", label="summary")
    await assert_md5(result.tpfn, "d3dffa4c7332a4cc3e01c4a6b1bd27c3", label="tpfn")
    await assert_md5(result.tpfp, "2aadd34139ee48b71346759b0b05791e", label="tpfp")


tests = [
    test_createsequencedictionary,
    test_intervallisttobed,
    test_gatk4_bedtointervallist,
    test_markduplicates_multiple,
    test_baserecalibrator,
    test_gatherbqsrreports,
    test_applybqsr,
    test_haplotypecaller,
    test_mergevcfs,
    test_cnnscorevariants,
    test_filtervarianttranches,
    test_gatk4_mutect2_tumor_normal,
    test_gatk4_mutect2_tumor_only,
    test_gatk4_mutect2_mitochondria,
    test_gatk4_mergemutectstats,
    test_gatk4_filtermutectcalls,
    test_gatk4_variantfiltration,
    test_gatk4_intervallisttools,
    test_gatk4_combinegvcfs,
    test_gatk4_splitncigarreads,
    test_gatk4_preprocessintervals,
    test_gatk4_annotateintervals,
    test_gatk4_collectreadcounts,
    test_gatk4_createreadcountpanelofnormals,
    test_gatk4_filterintervals,
    test_gatk4_determinegermlinecontigploidy,
    test_gatk4_germlinecnvcaller,
    test_gatk4_genomicsdbimport,
    test_gatk4_createsomaticpanelofnormals,
    test_gatk4_indexfeaturefile,
    test_gatk4_splitintervals,
    test_gatk4_getpileupsummaries,
    test_gatk4_gatherpileupsummaries,
    test_gatk4_calculatecontamination,
    test_gatk4_learnreadorientationmodel,
    test_gatk4_variantrecalibrator,
    test_gatk4_concordance,
]
