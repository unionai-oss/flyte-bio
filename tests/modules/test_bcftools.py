"""Tests for flyte_bio.modules.bcftools.

Upstream's bcftools/stats snapshot records the report's first six lines,
which include the command line bcftools records.

The other cases and expected values are upstream's bcftools/sort, norm, view,
filter, annotate, mpileup, merge, index, concat, reheader, consensus,
pluginsetgt, query, convert and pluginsplit tests
(modules commit 9a1cc9f, bcftools 1.25). Their snapshots hash outputs three ways, all reproduced here:

- a plain md5 (``assert_md5``), or for ``.gz`` outputs the md5 of the
  decompressed file (``assert_gunzipped_md5``). These cases pass
  ``--no-version``, so bcftools records no command line or date in the header;
- nft-vcf's ``variantsMD5`` (:func:`variants_md5`): the md5 of the record
  lines, concatenated without newlines;
- the file name only, for BCF outputs and VCFs whose header records a
  command line. Those get content checks here too.

Index files are snapshotted by name only.
"""

import asyncio
import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bcftools import (
    BcftoolsResult,
    bcftools_annotate,
    bcftools_concat,
    bcftools_consensus,
    bcftools_convert,
    bcftools_filter,
    bcftools_index,
    bcftools_merge,
    bcftools_mpileup,
    bcftools_norm,
    bcftools_pluginsetgt,
    bcftools_pluginsplit,
    bcftools_query,
    bcftools_reheader,
    bcftools_sort,
    bcftools_stats,
    bcftools_view,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

SC2 = "genomics/sarscov2/"
VCF = SC2 + "illumina/vcf/"
HS = "genomics/homo_sapiens/"
# test.vcf, which sort, view and filter (with --no-version) reproduce exactly.
TEST_VCF_MD5 = "8e722884ffb75155212a3fc053918766"
# Upstream's index configs: (label, args, expected index suffix).
INDEX_CASES = [
    ("vcf_gz_index", "--output-type z --write-index", "csi"),
    ("vcf_gz_index_csi", "--output-type z --write-index=csi", "csi"),
    ("vcf_gz_index_tbi", "--output-type z --write-index=tbi", "tbi"),
]


async def write_text(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


async def text_lines(file: File) -> list[str]:
    """A VCF's (or any text file's) lines, decompressing it if it's gzipped."""
    async with file.open("rb") as fh:
        raw = bytes(await fh.read())
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode().splitlines()


async def records(file: File) -> list[list[str]]:
    """A VCF's record lines, split into fields."""
    return [line.split("\t") for line in await text_lines(file) if not line.startswith("#")]


async def variants_md5(file: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of the record lines joined without newlines, in hex without leading zeros."""
    body = "".join(line for line in await text_lines(file) if not line.startswith("#"))
    return format(int(hashlib.md5(body.encode()).hexdigest(), 16), "x")


async def check(
    result: BcftoolsResult,
    name: str,
    label: str,
    md5: str | None = None,
    gunzipped_md5: str | None = None,
    variants: str | None = None,
    index: str | None = None,
) -> None:
    """Check an output's name, its index's name (None: no index) and whichever hash upstream snapshots."""
    assert result.vcf.name == name, f"{label}: output is {result.vcf.name}, expected {name}"
    actual_index = result.index.name if result.index else None
    assert actual_index == index, f"{label}: index is {actual_index}, expected {index}"
    if md5 is not None:
        await assert_md5(result.vcf, md5, label=label)
    if gunzipped_md5 is not None:
        await assert_gunzipped_md5(result.vcf, gunzipped_md5, label=label)
    if variants is not None:
        actual = await variants_md5(result.vcf)
        assert actual == variants, f"{label}: variantsMD5 {actual}, expected {variants}"


async def as_vcf(bcf: File) -> File:
    """A BCF converted to plain VCF, for content checks."""
    return (await bcftools_view(bcf, "converted", args="--output-type v")).vcf


@env.task
async def test_bcftools_stats() -> None:
    # upstream case: bcftools/stats "sarscov2 - vcf_gz"
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz")
    stats = await bcftools_stats(vcf, prefix="test")
    assert stats.path.endswith("/test.bcftools_stats.txt"), stats.path
    async with stats.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert lines[:6] == [
        "# This file was produced by bcftools stats (1.23.1+htslib-1.23.1) and can be plotted using plot-vcfstats.",
        "# The command line was:\tbcftools stats  test.vcf.gz",
        "#",
        "# Definition of sets:",
        "# ID\t[2]id\t[3]tab-separated file names",
        "ID\t0\ttest.vcf.gz",
    ], lines[:6]


@env.task
async def test_bcftools_sort_vcf() -> None:
    # upstream cases: bcftools/sort "sarscov2 - vcf", "... - vcf_gz_index", "... - vcf_gz_index_csi"
    # and "... - vcf_gz_index_tbi" (bcftools sort records no command line)
    vcf = await fixture(VCF + "test.vcf")
    plain, *indexed = await asyncio.gather(
        bcftools_sort(vcf, "test"), *(bcftools_sort(vcf, "test_vcf", args=args) for _, args, _ in INDEX_CASES)
    )
    await check(plain, "test.vcf.gz", "vcf", gunzipped_md5=TEST_VCF_MD5)
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(result, "test_vcf.vcf.gz", label, gunzipped_md5=TEST_VCF_MD5, index=f"test_vcf.vcf.gz.{suffix}")


@env.task
async def test_bcftools_norm_vcf() -> None:
    # upstream cases: bcftools/norm "sarscov2 - [ vcf, [] ], fasta" (and its three index variants),
    # "sarscov2 - [ vcf, tbi ], fasta", "... - vcf output", "... - vcf_gz output" (the same config as
    # "vcf output"), "... - bcf output" and "... - bcf_gz output"
    vcf, tbi, fasta = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"), fixture(VCF + "test.vcf.gz.tbi"), fixture(SC2 + "genome/genome.fasta.gz")
    )
    split = "-m -any --no-version"
    no_index, with_index, vcf_out, bcf_gz_out, bcf_out, *indexed = await asyncio.gather(
        bcftools_norm(vcf, fasta, "test_norm", args=split),
        bcftools_norm(vcf, fasta, "test_norm", index=tbi, args=split),
        bcftools_norm(vcf, fasta, "test_norm", index=tbi, args=f"{split} --output-type v"),
        bcftools_norm(vcf, fasta, "test_norm", index=tbi, args=f"{split} --output-type b"),
        bcftools_norm(vcf, fasta, "test_norm", index=tbi, args=f"{split} --output-type u"),
        *(bcftools_norm(vcf, fasta, "test_vcf", args=f"{args} --no-version") for _, args, _ in INDEX_CASES),
    )
    normalized = "63e5adbaf3dd94550e9e3d7935dd28db"
    await check(no_index, "test_norm.vcf.gz", "[vcf, []]", gunzipped_md5=normalized)
    await check(with_index, "test_norm.vcf.gz", "[vcf, tbi]", gunzipped_md5=normalized)
    await check(vcf_out, "test_norm.vcf", "vcf output", md5=normalized)
    await check(bcf_gz_out, "test_norm.bcf.gz", "bcf output", gunzipped_md5="638c3c25bdd495c90ecbccb69ee77f07")
    await check(bcf_out, "test_norm.bcf", "bcf_gz output", md5="bf88706ef69c44ca9e287bc953ba3593")
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(result, "test_vcf.vcf.gz", label, gunzipped_md5=normalized, index=f"test_vcf.vcf.gz.{suffix}")


@env.task
async def test_bcftools_view_vcf() -> None:
    # upstream cases: bcftools/view "sarscov2 - [vcf, tbi], [], [], []" (and its three index variants)
    # and "sarscov2 - [vcf, tbi], vcf, tsv, []"
    vcf, tbi, regions, targets = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"),
        fixture(VCF + "test.vcf.gz.tbi"),
        fixture(VCF + "test3.vcf.gz"),
        fixture(VCF + "test2.targets.tsv.gz"),
    )
    args = "--no-version --output-type v"
    plain, subset, *indexed = await asyncio.gather(
        bcftools_view(vcf, "out", index=tbi, args=args),
        bcftools_view(vcf, "out", index=tbi, regions=regions, targets=targets, args=args),
        *(bcftools_view(vcf, "out_vcf", index=tbi, args=f"{a} --no-version") for _, a, _ in INDEX_CASES),
    )
    await check(plain, "out.vcf", "[vcf, tbi]", md5=TEST_VCF_MD5)
    await check(subset, "out.vcf", "regions and targets", md5="1bcbd0eff25d316ba915d06463aab17b")
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(result, "out_vcf.vcf.gz", label, gunzipped_md5=TEST_VCF_MD5, index=f"out_vcf.vcf.gz.{suffix}")


@env.task
async def test_bcftools_filter_expression() -> None:
    # Not an upstream case: an expression with spaces reaches bcftools whole, via include/exclude.
    vcf = await fixture(VCF + "test.vcf")
    every = await records(vcf)
    cutoff = sorted(int(r[1]) for r in every)[len(every) // 2]
    kept = [r for r in every if int(r[1]) > cutoff and r[3] != "N"]
    assert 0 < len(kept) < len(every), "the test VCF should be split by the expression"

    filtered, viewed = await asyncio.gather(
        bcftools_filter(vcf, "expr_filter", include=f"POS > {cutoff} && REF != 'N'", args="--no-version -Ov"),
        bcftools_view(vcf, "expr_view", exclude=f"POS <= {cutoff} || REF == 'N'", args="--no-version -Ov"),
    )
    assert [r[1] for r in await records(filtered.vcf)] == [r[1] for r in kept], "filter --include"
    assert [r[1] for r in await records(viewed.vcf)] == [r[1] for r in kept], "view --exclude"


@env.task
async def test_bcftools_filter_vcf() -> None:
    # upstream cases: bcftools/filter "sarscov2 - vcf" (and its three index variants),
    # "sarscov2 - vcf - bcf output" and "sarscov2 - vcf.gz, tbi - region filter"
    vcf, vcf_gz, tbi = await asyncio.gather(
        fixture(VCF + "test.vcf"), fixture(VCF + "test.vcf.gz"), fixture(VCF + "test.vcf.gz.tbi")
    )
    plain, bcf, region, *indexed = await asyncio.gather(
        bcftools_filter(vcf, "vcf_test", args="--no-version --output-type z"),
        bcftools_filter(vcf, "bcf_test", args="--no-version --output-type b"),
        bcftools_filter(vcf_gz, "bcf_test_vcf", index=tbi, args="--output-type z -r MT192765.1 --no-version"),
        *(bcftools_filter(vcf, "vcf_test_vcf", args=f"{args} --no-version") for _, args, _ in INDEX_CASES),
    )
    await check(plain, "vcf_test.vcf.gz", "vcf", gunzipped_md5=TEST_VCF_MD5)
    await check(bcf, "bcf_test.bcf.gz", "bcf output", gunzipped_md5="c8a304c8d2892039201154153c8cd536")
    await check(region, "bcf_test_vcf.vcf.gz", "region filter", gunzipped_md5=TEST_VCF_MD5)
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(
            result, "vcf_test_vcf.vcf.gz", label, gunzipped_md5=TEST_VCF_MD5, index=f"vcf_test_vcf.vcf.gz.{suffix}"
        )


@env.task
async def test_bcftools_annotate_vcf() -> None:
    # upstream cases: bcftools/annotate "sarscov2 - [vcf, tbi, annotation, annotation_tbi], [], [], [] -
    # vcf_output", the same without tbi, and the three index variants
    vcf, tbi, annotations, annotations_tbi = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"),
        fixture(VCF + "test.vcf.gz.tbi"),
        fixture(VCF + "test2.vcf.gz"),
        fixture(VCF + "test2.vcf.gz.tbi"),
    )
    strip = "-x ID,INFO/DP,FORMAT/DP --output-type z"
    with_index, no_index, *indexed = await asyncio.gather(
        bcftools_annotate(vcf, "test_ann", tbi, annotations, annotations_tbi, args=strip),
        bcftools_annotate(vcf, "test_ann", None, annotations, annotations_tbi, args=strip),
        *(
            bcftools_annotate(vcf, "test_ann", tbi, annotations, annotations_tbi, args=f"{args} --no-version")
            for _, args, _ in INDEX_CASES
        ),
    )
    stripped = "250b64289ab9d48f76359d01699fdf7d"
    await check(with_index, "test_ann.vcf.gz", "vcf_output", variants=stripped)
    await check(no_index, "test_ann.vcf.gz", "vcf_output without tbi", variants=stripped)
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(
            result,
            "test_ann.vcf.gz",
            label,
            variants="bc7bf3ee9e8430e064c539eb81e59bf9",
            index=f"test_ann.vcf.gz.{suffix}",
        )


@env.task
async def test_bcftools_annotate_files() -> None:
    # upstream cases: bcftools/annotate "sarscov2 - [vcf, [], annotation, annotation_tbi], [], header, [] -
    # bcf_output", "... columns, [], [] - bcf_output" and "... [], header, rename_chrs - vcf_gz_index".
    # Their snapshots record only names; the content checks are ours.
    headers = (
        '##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">\n'
        '##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">\n'
    )
    vcf, annotations, annotations_tbi, header_lines, columns, rename = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"),
        fixture(VCF + "test2.vcf.gz"),
        fixture(VCF + "test2.vcf.gz.tbi"),
        write_text("headers.vcf", headers),
        write_text("columns.txt", "INFO/ICB\nINFO/HOB\nINFO/DP4\n"),
        write_text("rename.txt", "MT192765.1 renamed\n"),
    )
    strip = "-x ID,INFO/DP,FORMAT/DP --output-type u"
    with_header, with_columns, renamed = await asyncio.gather(
        bcftools_annotate(vcf, "test_ann", None, annotations, annotations_tbi, header_lines=header_lines, args=strip),
        bcftools_annotate(vcf, "test_ann", None, annotations, annotations_tbi, columns=columns, args=strip),
        bcftools_annotate(
            vcf,
            "test_ann",
            None,
            annotations,
            annotations_tbi,
            header_lines=header_lines,
            rename_chrs=rename,
            args="--output-type z --write-index --no-version",
        ),
    )
    await check(with_header, "test_ann.bcf", "header bcf_output")
    await check(with_columns, "test_ann.bcf", "columns bcf_output")
    await check(renamed, "test_ann.vcf.gz", "rename_chrs vcf_gz_index", index="test_ann.vcf.gz.csi")

    original = await records(vcf)
    for label, result in (("header", with_header), ("columns", with_columns)):
        lines = await text_lines(await as_vcf(result.vcf))
        body = [line.split("\t") for line in lines if not line.startswith("#")]
        assert len(body) == len(original), f"{label}: {len(body)} records, expected {len(original)}"
        for record in body:
            info_keys = [field.split("=")[0] for field in record[7].split(";")]
            assert record[2] == "." and "DP" not in info_keys, f"{label}: ID or INFO/DP kept: {record[:8]}"
        if label == "header":
            assert all(h in lines for h in headers.splitlines()), f"{label}: header lines missing"
        else:
            indel = next(r for r in body if r[1] == "10506")
            assert "ICB=1" in indel[7] and "HOB=0.5" in indel[7], f"columns: ICB/HOB not annotated: {indel[7]}"
    renamed_records = await records(renamed.vcf)
    assert len(renamed_records) == len(original), f"rename_chrs: {len(renamed_records)} records"
    assert {r[0] for r in renamed_records} == {"renamed"}, {r[0] for r in renamed_records}


async def query_targets() -> File:
    """Upstream's setup step, ``bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n'`` of test.vcf, as ``test.txt``."""
    rows = await records(await fixture(VCF + "test.vcf"))
    return await write_text("test.txt", "".join(f"{r[0]}\t{r[1]}\t{r[3]},{r[4]}\n" for r in rows))


@env.task
async def test_bcftools_mpileup_bam() -> None:
    # upstream cases: bcftools/mpileup "sarscov2 - [bam, []], fasta, false", "... true",
    # "[bam, bed_mpileup], fasta, false", "[bam, bed_call], fasta, false" and
    # "[bam, bed_mpileup, bed_call], fasta, false"; the targets file comes from upstream's bcftools/query
    # setup step, recomputed by query_targets
    bam, fasta, bed, targets = await asyncio.gather(
        fixture(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(SC2 + "genome/genome.fasta"),
        fixture(SC2 + "genome/bed/test.bed"),
        query_targets(),
    )
    args2 = "--no-version --ploidy 1 --multiallelic-caller"
    plain, saved, bed_mpileup, bed_call, both = await asyncio.gather(
        bcftools_mpileup(bam, fasta, "test", args2=args2, args3="--no-version --write-index=csi --output-type z"),
        bcftools_mpileup(
            bam, fasta, "test", save_mpileup=True, args2=args2, args3="--no-version --write-index=tbi --output-type z"
        ),
        bcftools_mpileup(bam, fasta, "test", intervals_mpileup=bed, args2=args2),
        bcftools_mpileup(bam, fasta, "test", intervals_call=targets, args2=args2),
        bcftools_mpileup(bam, fasta, "test", intervals_mpileup=targets, intervals_call=targets, args2=args2),
    )
    cases = [
        ("[bam, []] false", plain, "test.vcf.gz", "test.vcf.gz.csi", "113f3e0b28ca8d7746963bd5a0d40bc0",
         "e534df91b8c27be0eb520bd02182d62b"),
        ("[bam, []] true", saved, "test.vcf.gz", "test.vcf.gz.tbi", "113f3e0b28ca8d7746963bd5a0d40bc0",
         "e534df91b8c27be0eb520bd02182d62b"),
        ("bed_mpileup", bed_mpileup, "test.vcf", None, "7341fd584a827cddfc5f181a81ed5076",
         "2d94085e36c72c784c1d077d41a8bd4e"),
        ("bed_call", bed_call, "test.vcf", None, "46cbd80967ad92b1d603492e503e6a06",
         "5faeeeedc7308bd2d8596ec9b864dafa"),
        ("bed_mpileup, bed_call", both, "test.vcf", None, "c3e829c86ad6ce72945789383bb0108f",
         "76808a644b9d73ee77c15795e61f679a"),
    ]  # fmt: skip
    for label, result, name, index, variants, stats_md5 in cases:
        assert result.vcf.name == name, f"{label}: output is {result.vcf.name}, expected {name}"
        actual_index = result.index.name if result.index else None
        assert actual_index == index, f"{label}: index is {actual_index}, expected {index}"
        actual = await variants_md5(result.vcf)
        assert actual == variants, f"{label}: variantsMD5 {actual}, expected {variants}"
        await assert_md5(result.stats, stats_md5, label=f"{label} stats")
        assert (result.mpileup is not None) == (result is saved), f"{label}: mpileup is {result.mpileup}"
    assert saved.mpileup is not None and saved.mpileup.name == "test.mpileup.gz", saved.mpileup
    await assert_gunzipped_md5(saved.mpileup, "4dd120039fcf4ff945d0ca749ad9465b", label="mpileup")


@env.task
async def test_bcftools_merge_vcfs() -> None:
    # upstream cases: bcftools/merge "sarscov2 - [vcf, tbi], [], [], []" with the default, vcf, vcf.gz,
    # bcf, bcf.gz and the three index outputs, "... [], [], bed" and "... - one sample"
    vcfs = list(await asyncio.gather(fixture(VCF + "test2.vcf.gz"), fixture(VCF + "test3.vcf.gz")))
    tbis = list(await asyncio.gather(fixture(VCF + "test2.vcf.gz.tbi"), fixture(VCF + "test3.vcf.gz.tbi")))
    bed = await fixture(SC2 + "genome/bed/test.bed")
    forced = "--force-samples --force-single --no-version"
    default, vcf_out, vcf_gz, bcf, bcf_gz, regions, one_sample, *indexed = await asyncio.gather(
        bcftools_merge(vcfs, tbis, "test", args=forced),
        bcftools_merge(vcfs, tbis, "test", args="--output-type v --no-version"),
        bcftools_merge(vcfs, tbis, "test", args="--output-type z --no-version"),
        bcftools_merge(vcfs, tbis, "test", args="--output-type u --no-version"),
        bcftools_merge(vcfs, tbis, "test", args="--output-type b --no-version"),
        bcftools_merge(vcfs, tbis, "test", bed=bed, args=forced),
        bcftools_merge(vcfs[:1], tbis[:1], "test", args=forced),
        *(bcftools_merge(vcfs, tbis, "test", args=f"{args} --no-version") for _, args, _ in INDEX_CASES),
    )
    merged = "e0de448dc8e712956a03ce68d79a0b3a"
    await check(default, "test.vcf", "default", variants=merged)
    await check(vcf_out, "test.vcf", "vcf output", md5="57bb84274f336465d0a0946b532093b0")
    await check(vcf_gz, "test.vcf.gz", "vcf.gz output", variants=merged)
    await check(bcf, "test.bcf", "bcf output")
    await check(bcf_gz, "test.bcf.gz", "bcf.gz output")
    await check(regions, "test.vcf", "bed", md5="febdcfb851dcfc83d8248520830aef10")
    await check(one_sample, "test.vcf", "one sample", md5="2a374cf02f0c32cf607646167e7f153b")
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(result, "test.vcf.gz", label, variants=merged, index=f"test.vcf.gz.{suffix}")
    # The BCF snapshots record only names; the merged records must match the VCF's.
    for label, result in (("bcf", bcf), ("bcf.gz", bcf_gz)):
        actual = await variants_md5(await as_vcf(result.vcf))
        assert actual == merged, f"{label} output: variantsMD5 {actual}, expected {merged}"


@env.task
async def test_bcftools_merge_gvcf() -> None:
    # upstream case: bcftools/merge "homo_sapiens - [vcf, tbi], fasta, fai, bed  - vcf.gz output"
    gvcf = HS + "illumina/gvcf/"
    vcfs, tbis, bed, fasta, fai = await asyncio.gather(
        asyncio.gather(fixture(gvcf + "test.genome.vcf.gz"), fixture(gvcf + "test2.genome.vcf.gz")),
        asyncio.gather(fixture(gvcf + "test.genome.vcf.gz.tbi"), fixture(gvcf + "test2.genome.vcf.gz.tbi")),
        fixture(HS + "genome/genome.bed"),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/genome.fasta.fai"),
    )
    result = await bcftools_merge(
        list(vcfs),
        list(tbis),
        "test",
        bed=bed,
        fasta=fasta,
        fai=fai,
        args="--force-samples --no-version --output-type z --gvcf genome.fasta",
    )
    await check(result, "test.vcf.gz", "gvcf", variants="645b7f7f9131bfe350a9ec3cf82c17fe")


@env.task
async def test_bcftools_index_vcf() -> None:
    # upstream cases: bcftools/index "sarscov2 - vcf - csi" and "sarscov2 - vcf - tbi" (names only);
    # that each index serves a region query is our check
    vcf = await fixture(VCF + "test.vcf.gz")
    csi, tbi = await asyncio.gather(bcftools_index(vcf), bcftools_index(vcf, args="--tbi"))
    for index, name in ((csi, "test.vcf.gz.csi"), (tbi, "test.vcf.gz.tbi")):
        assert index is not None and index.name == name, f"index is {index and index.name}, expected {name}"
    regions = await asyncio.gather(
        *(
            bcftools_filter(vcf, f"region_{suffix}", index=index, args="--output-type z -r MT192765.1 --no-version")
            for suffix, index in (("csi", csi), ("tbi", tbi))
        )
    )
    for suffix, result in zip(("csi", "tbi"), regions, strict=True):
        await check(result, f"region_{suffix}.vcf.gz", f"region query with the {suffix}", gunzipped_md5=TEST_VCF_MD5)


@env.task
async def test_bcftools_concat_vcfs() -> None:
    # upstream cases: bcftools/concat "homo_sapiens - [[vcf1, vcf2], [tbi1, tbi2]]" (and its three index
    # variants) and "homo_sapiens - [[vcf1, vcf2], []]"
    hs_vcf = HS + "illumina/vcf/"
    vcfs = list(
        await asyncio.gather(
            fixture(hs_vcf + "NA12878.chr21_22.1X.glimpse2.vcf.gz"), fixture(hs_vcf + "NA12878_GIAB.chr21_22.vcf.gz")
        )
    )
    indexes = list(
        await asyncio.gather(
            fixture(hs_vcf + "NA12878_GIAB.chr21_22.vcf.gz.csi"),
            fixture(hs_vcf + "NA12878.chr21_22.1X.glimpse2.vcf.gz.csi"),
        )
    )
    index_cases = [
        ("vcf_gz_index", "--output-type z --write-index --no-version --allow-overlaps", "csi"),
        ("vcf_gz_index_csi", "--output-type z --write-index=csi --no-version --allow-overlaps -a", "csi"),
        ("vcf_gz_index_tbi", "--output-type z --write-index=tbi --no-version --allow-overlaps", "tbi"),
    ]
    plain, unindexed, *indexed = await asyncio.gather(
        bcftools_concat(vcfs, indexes, "test3", args="--no-version -a"),
        bcftools_concat(vcfs, [], "test3", args="--no-version -a"),
        *(bcftools_concat(vcfs, indexes, "test3_vcf", args=args) for _, args, _ in index_cases),
    )
    concatenated = "0b4986bbb6bdc9e0ee12ba09e70092e9"
    await check(plain, "test3.vcf", "[vcfs, indexes]", md5=concatenated)
    await check(unindexed, "test3.vcf", "[vcfs, []]", md5=concatenated)
    for (label, _, suffix), result in zip(index_cases, indexed, strict=True):
        await check(result, "test3_vcf.vcf.gz", label, gunzipped_md5=concatenated, index=f"test3_vcf.vcf.gz.{suffix}")


@env.task
async def test_bcftools_reheader_vcf() -> None:
    # upstream cases: bcftools/reheader "sarscov2 - [vcf, [], []], fai - vcf output", "... - vcf.gz output",
    # "... - bcf output", the three index variants, "sarscov2 - [vcf, header, []], []" and
    # "sarscov2 - [vcf, [], samples], fai"
    vcf, fai, header, samples = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"),
        fixture(SC2 + "genome/genome.fasta.fai"),
        fixture(VCF + "test.vcf"),
        write_text("samples.txt", "new_name\n"),
    )
    vcf_out, vcf_gz, bcf, with_header, renamed, *indexed = await asyncio.gather(
        bcftools_reheader(vcf, "tested", fai=fai, args2="--no-version"),
        bcftools_reheader(vcf, "tested", fai=fai, args2="--no-version --output-type z"),
        bcftools_reheader(vcf, "tested", fai=fai, args2="--no-version --output-type b"),
        bcftools_reheader(vcf, "tested", header=header, args2="--no-version"),
        bcftools_reheader(vcf, "tested", samples=samples, fai=fai, args2="--no-version"),
        *(bcftools_reheader(vcf, "test_vcf", fai=fai, args2=f"{args} --no-version") for _, args, _ in INDEX_CASES),
    )
    await check(vcf_out, "tested.vcf", "vcf output", md5=TEST_VCF_MD5)
    await check(vcf_gz, "tested.vcf.gz", "vcf.gz output", gunzipped_md5=TEST_VCF_MD5)
    await check(bcf, "tested.bcf.gz", "bcf output", gunzipped_md5="c8a304c8d2892039201154153c8cd536")
    await check(with_header, "tested.vcf", "header", md5="3189bc9a720d5d5d3006bf72d91300cb")
    await check(renamed, "tested.vcf", "samples", md5="c64c373c10b0be24b29d6f18708ec1e8")
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(result, "test_vcf.vcf.gz", label, gunzipped_md5=TEST_VCF_MD5, index=f"test_vcf.vcf.gz.{suffix}")


@env.task
async def test_bcftools_consensus_fasta() -> None:
    # upstream cases: bcftools/consensus "bcftools - test" (masked by test.bed) and "bcftools - test - no mask"
    vcf, tbi, fasta, bed = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"),
        fixture(VCF + "test.vcf.gz.tbi"),
        fixture(SC2 + "genome/genome.fasta"),
        fixture(SC2 + "genome/bed/test.bed"),
    )
    masked, unmasked = await asyncio.gather(
        bcftools_consensus(vcf, tbi, fasta, "test", mask=bed), bcftools_consensus(vcf, tbi, fasta, "test")
    )
    for label, consensus, md5 in (
        ("masked", masked, "313f4a657187fdaf293c4ed69d98c112"),
        ("no mask", unmasked, "e57d3e4d113f989bc069c5fd61627091"),
    ):
        assert consensus.name == "test.fa", f"{label}: output is {consensus.name}"
        await assert_md5(consensus, md5, label=f"consensus {label}")


@env.task
async def test_bcftools_pluginsetgt_vcf() -> None:
    # upstream cases: bcftools/pluginsetgt "sarscov2 - [vcf, tbi], [], []" (and its three index variants)
    # and "sarscov2 - [vcf, tbi], vcf, tsv, []": every genotype phased, with the plugin's -s 42
    vcf, tbi, regions, targets = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"),
        fixture(VCF + "test.vcf.gz.tbi"),
        fixture(VCF + "test3.vcf.gz"),
        fixture(VCF + "test2.targets.tsv.gz"),
    )
    common = {"target_gt": "a", "new_gt": "p", "index": tbi, "args2": "-s 42"}
    args = "--no-version --output-type v"
    plain, subset, *indexed = await asyncio.gather(
        bcftools_pluginsetgt(vcf, "out", args=args, **common),
        bcftools_pluginsetgt(vcf, "out", regions=regions, targets=targets, args=args, **common),
        *(bcftools_pluginsetgt(vcf, "out_vcf", args=f"{a} --no-version", **common) for _, a, _ in INDEX_CASES),
    )
    phased = "9fb75c9fdb304b49c3895476666a85bb"
    await check(plain, "out.vcf", "[vcf, tbi]", md5=phased)
    await check(subset, "out.vcf", "regions and targets", md5="55e955c359bf9094bb692660f3c1c677")
    for (label, _, suffix), result in zip(INDEX_CASES, indexed, strict=True):
        await check(result, "out_vcf.vcf.gz", label, gunzipped_md5=phased, index=f"out_vcf.vcf.gz.{suffix}")


@env.task
async def test_bcftools_query_vcf() -> None:
    # upstream cases: bcftools/query "sarscov2 - [vcf, tbi], [], [], []" and "... vcf, tsv, []",
    # with the test config's -f '%CHROM %POS %REF %ALT[%SAMPLE=%GT]'
    vcf, tbi, regions, targets = await asyncio.gather(
        fixture(VCF + "test.vcf.gz"),
        fixture(VCF + "test.vcf.gz.tbi"),
        fixture(VCF + "test3.vcf.gz"),
        fixture(VCF + "test2.targets.tsv.gz"),
    )
    query_format = "%CHROM %POS %REF %ALT[%SAMPLE=%GT]"
    plain, subset = await asyncio.gather(
        bcftools_query(vcf, "out", index=tbi, format=query_format),
        bcftools_query(vcf, "out", index=tbi, regions=regions, targets=targets, format=query_format),
    )
    for label, output, md5 in (
        ("[vcf, tbi]", plain, "87a2ab194e1ee3219b44e58429ec3307"),
        ("regions and targets", subset, "75a6bd0084e2e1838cf7baba11b99d19"),
    ):
        assert output.name == "out.txt", f"{label}: output is {output.name}"
        await assert_md5(output, md5, label=f"query {label}")


@env.task
async def test_bcftools_convert_gvcf() -> None:
    # upstream cases: bcftools/convert "homo_sapiens - vcf1, tbi, [], fasta - to bcf", "... - to hap"
    # and "homo_sapiens - vcf1, tbi, bed, fasta - to gvcf"
    gvcf, tbi, fasta, bed = await asyncio.gather(
        fixture(HS + "illumina/gvcf/test.genome.vcf.gz"),
        fixture(HS + "illumina/gvcf/test.genome.vcf.gz.tbi"),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/genome.bed"),
    )
    to_bcf, to_hap, to_vcf = await asyncio.gather(
        bcftools_convert(
            gvcf, "test3", index=tbi, fasta=fasta, args="--no-version --gvcf2vcf --output-type b5 --write-index=tbi"
        ),
        bcftools_convert(gvcf, "test3", index=tbi, fasta=fasta, args="--no-version --haplegendsample test"),
        bcftools_convert(
            gvcf, "test3", index=tbi, fasta=fasta, bed=bed, args="--no-version --gvcf2vcf --write-index=csi"
        ),
    )
    # bcftools indexes a BCF with a .csi whatever --write-index asks for, as upstream's snapshot shows.
    for label, result, name, md5 in (
        ("to bcf", to_bcf, "test3.bcf.gz", "3e26eb6a474996a4ecd6fb874ff4dc97"),
        ("to gvcf", to_vcf, "test3.vcf.gz", "ef40139c7f15d8df22c0cc5cb409ef88"),
    ):
        assert result.vcf is not None and result.vcf.name == name, f"{label}: output is {result.vcf}"
        assert result.index is not None and result.index.name == f"{name}.csi", f"{label}: index is {result.index}"
        await assert_gunzipped_md5(result.vcf, md5, label=f"convert {label}")
    assert to_hap.hap is not None and to_hap.legend is not None and to_hap.samples is not None
    await assert_gunzipped_md5(to_hap.hap, "1ecdaca6b10a55409ffd245fea63014f", label="test.hap.gz")
    await assert_gunzipped_md5(to_hap.legend, "2b537fa86c2d86b175ebbd2007e05371", label="test.legend.gz")
    await assert_md5(to_hap.samples, "2606935b2fe508d5fd14afc249c141c0", label="test.samples")


SPLIT_VCF = HS + "illumina/gatk/paired_mutect2_calls/test_test2_paired_mutect2_calls.vcf.gz"
# Upstream's samples and groups files: each sample kept under its own name, in its own file.
SPLIT_SAMPLES = "normal\t-\tnormal\ntumour\t-\ttumour\n"


@env.task
async def test_bcftools_pluginsplit_vcf() -> None:
    # upstream cases: bcftools/pluginsplit "homo_sapiens - [ vcf, tbi ], samples, [], [], []",
    # "... [], groups, regions, targets" and its "- tbi" variant (with the test config's args)
    vcf, tbi, regions, targets, samples, groups = await asyncio.gather(
        fixture(SPLIT_VCF),
        fixture(SPLIT_VCF + ".tbi"),
        fixture(HS + "genome/genome.bed"),
        fixture(HS + "genome/genome.multi_intervals.bed"),
        write_text("samples.txt", SPLIT_SAMPLES),
        write_text("groups.txt", SPLIT_SAMPLES),
    )
    by_samples, by_groups, indexed = await asyncio.gather(
        bcftools_pluginsplit(vcf, index=tbi, samples=samples),
        bcftools_pluginsplit(vcf, index=tbi, groups=groups, regions=regions, targets=targets),
        bcftools_pluginsplit(
            vcf, index=tbi, groups=groups, regions=regions, targets=targets, args="--write-index=tbi --output-type z"
        ),
    )
    for label, result, md5s in (
        ("samples", by_samples, ["69cfc4bf92bf3e2847081a2026a4d3bb", "08fa5c8d5561c2a8d7c300cb0eea1042"]),
        ("groups", by_groups, ["78e962c3a724790b27041f28a4180600", "3008a0a06a2f319897fa04ea0018fcdd"]),
    ):
        names = [f.name for f in result.vcfs]
        assert names == ["normal.vcf", "tumour.vcf"], f"{label}: wrote {names}"
        for output, md5 in zip(result.vcfs, md5s, strict=True):
            await assert_md5(output, md5, label=f"split {label} {output.name}")
    names = [f.name for f in indexed.vcfs + indexed.indexes]
    assert names == ["normal.vcf.gz", "tumour.vcf.gz", "normal.vcf.gz.tbi", "tumour.vcf.gz.tbi"], names


@env.task
async def test_bcftools_pluginsplit_no_samples() -> None:
    # upstream case: bcftools/pluginsplit "homo_sapiens - [ vcf, tbi ], [], [], [], [], - error no sample"
    vcf = await fixture(HS + "genome/vcf/dbsnp_146.hg38.vcf.gz")
    tbi = await fixture(HS + "genome/vcf/dbsnp_146.hg38.vcf.gz.tbi")
    try:
        await bcftools_pluginsplit(vcf, index=tbi, args="--write-index=tbi --output-type z")
    except Exception as e:  # the task fails: the VCF has no samples
        # Upstream also matches bcftools' "No samples to split: input/dbsnp_146.hg38.vcf.gz"; a shell
        # task's stderr isn't part of the error the cluster reports, so only the failure is checked.
        print(f"bcftools plugin split failed as expected: {e}")
    else:
        raise AssertionError("bcftools plugin split accepted a VCF without samples")


tests = [
    test_bcftools_stats,
    test_bcftools_filter_expression,
    test_bcftools_sort_vcf,
    test_bcftools_norm_vcf,
    test_bcftools_view_vcf,
    test_bcftools_filter_vcf,
    test_bcftools_annotate_vcf,
    test_bcftools_annotate_files,
    test_bcftools_mpileup_bam,
    test_bcftools_merge_vcfs,
    test_bcftools_merge_gvcf,
    test_bcftools_index_vcf,
    test_bcftools_concat_vcfs,
    test_bcftools_reheader_vcf,
    test_bcftools_consensus_fasta,
    test_bcftools_pluginsetgt_vcf,
    test_bcftools_query_vcf,
    test_bcftools_convert_gvcf,
    test_bcftools_pluginsplit_vcf,
    test_bcftools_pluginsplit_no_samples,
]
