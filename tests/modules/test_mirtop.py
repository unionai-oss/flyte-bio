"""Tests for flyte_bio.modules.mirtop.

The cases and expected md5s are upstream's mirtop/gff, mirtop/counts,
mirtop/export and mirtop/stats tests (modules commit 9a1cc9f, mirtop 0.4.30),
with upstream's ``--add-extra`` for counts.
"""

from flyte.io import File

from flyte_bio.modules.mirtop import mirtop_counts, mirtop_export, mirtop_gff, mirtop_stats
from tests.framework import assert_md5, env, fixture

MIRTOP_DATA = "delete_me/mirtop/"
STATS_GFF = "https://github.com/miRTop/mirtop/raw/master/data/examples/gff/correct_file.gff"


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def references() -> tuple[File, File]:
    return await fixture(MIRTOP_DATA + "hairpin_mirtop.fa"), await fixture(MIRTOP_DATA + "hsa.gff3")


async def isomir_gff() -> File:
    hairpin, gtf = await references()
    bam = await fixture(MIRTOP_DATA + "sim_isomir_sort.bam")
    return await mirtop_gff([bam], hairpin, gtf, "hsa", prefix="sample_sim_isomir_bam")


@env.task
async def test_mirtop_gff_isomir_bam() -> None:
    gff = await isomir_gff()
    assert file_name(gff) == "sample_sim_isomir_bam_mirtop.gff", gff.path
    await assert_md5(gff, "da04f476c3fb8670e861fd8dd83418f9", label="sample_sim_isomir_bam_mirtop.gff")


@env.task
async def test_mirtop_counts_isomir_bam() -> None:
    hairpin, gtf = await references()
    tsv = await mirtop_counts(await isomir_gff(), hairpin, gtf, "hsa", args="--add-extra")
    assert file_name(tsv) == "sample_sim_isomir_bam_mirtop.tsv", tsv.path
    await assert_md5(tsv, "43f8a525104c2d9b5a8937564c3a14f6", label="sample_sim_isomir_bam_mirtop.tsv")


@env.task
async def test_mirtop_export_isomir_bam() -> None:
    hairpin, gtf = await references()
    result = await mirtop_export(await isomir_gff(), hairpin, gtf, "hsa")
    assert result.fasta == [] and result.vcf == [], "isomir export wrote fasta or vcf"
    names = [file_name(f) for f in result.tsv]
    assert names == ["sample_sim_isomir_bam_mirtop_rawData.tsv"], names
    await assert_md5(result.tsv[0], "efbcbe67716a4a56f89e538af2251dcc", label=names[0])


@env.task
async def test_mirtop_stats_correct_file() -> None:
    result = await mirtop_stats(await fixture(STATS_GFF), prefix="mirtop_gff_sample1")
    await assert_md5(result.txt, "e8e2307d0029bea45581765f6c469840", label="mirtop_gff_sample1_stats.txt")
    await assert_md5(result.log, "7325c65f141f53a891de72337f49c846", label="mirtop_gff_sample1_stats.log")


tests = [
    test_mirtop_gff_isomir_bam,
    test_mirtop_counts_isomir_bam,
    test_mirtop_export_isomir_bam,
    test_mirtop_stats_correct_file,
]
