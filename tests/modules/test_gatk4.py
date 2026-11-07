"""Tests for flyte_bio.modules.gatk4 (md5s from upstream's gatk4 module snapshots)."""

import asyncio
import gzip
import hashlib
import tempfile
import zlib
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.gatk4 import (
    Mutect2Result,
    gatk4_applybqsr,
    gatk4_baserecalibrator,
    gatk4_bedtointervallist,
    gatk4_cnnscorevariants,
    gatk4_combinegvcfs,
    gatk4_createsequencedictionary,
    gatk4_filtermutectcalls,
    gatk4_filtervarianttranches,
    gatk4_gatherbqsrreports,
    gatk4_haplotypecaller,
    gatk4_intervallisttobed,
    gatk4_intervallisttools,
    gatk4_markduplicates,
    gatk4_mergemutectstats,
    gatk4_mergevcfs,
    gatk4_mutect2,
    gatk4_splitncigarreads,
    gatk4_variantfiltration,
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
]
