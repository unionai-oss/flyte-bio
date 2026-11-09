"""Tests for flyte_bio.modules.popscle.

The cases are upstream's popscle/dscpileup, popscle/demuxlet and
popscle/freemuxlet tests (modules commit 9a1cc9f), with their test configs'
``ext.args`` (``--field GT`` for demuxlet, ``--aux-files`` for freemuxlet's
auxiliary case). Upstream snapshots only the outputs' names, so these tests
check each one exists and isn't empty.
"""

from flyte.io import Dir, File

from flyte_bio.modules.popscle import popscle_demuxlet, popscle_dscpileup, popscle_freemuxlet
from tests.framework import assert_nonempty, env, fixture

DATA = "genomics/homo_sapiens/demultiplexing/"


async def popscle_pileup() -> Dir:
    bam = await fixture(DATA + "chr21.bam")
    return await popscle_dscpileup(bam, await fixture(DATA + "donor_genotype_chr21.vcf"), prefix="sample1")


async def check_named(file: File | None, name: str) -> None:
    assert file is not None, f"no {name}"
    assert file.name == name, f"expected {name}, got {file.name}"
    await assert_nonempty(file, label=name)


@env.task
async def test_popscle_dscpileup() -> None:
    # upstream case: popscle/dscpileup "demultiplexing - bam"
    pileup = await popscle_pileup()
    for ext in ("cel", "plp", "var", "umi"):
        name = f"sample1.{ext}.gz"
        await check_named(await pileup.get_file(name), name)


@env.task
async def test_popscle_demuxlet_bam() -> None:
    # upstream case: popscle/demuxlet "demultiplexing - bam"
    best = await popscle_demuxlet(
        await fixture(DATA + "donor_genotype_chr21.vcf"),
        prefix="sample1",
        bam=await fixture(DATA + "chr21.bam"),
        args="--field GT",
    )
    await check_named(best, "sample1.best")


@env.task
async def test_popscle_demuxlet_pileup() -> None:
    # upstream case: popscle/demuxlet "demultiplexing - bam - pileup"
    best = await popscle_demuxlet(
        await fixture(DATA + "donor_genotype_chr21.vcf"),
        prefix="sample1",
        pileup=await popscle_pileup(),
        args="--field GT",
    )
    await check_named(best, "sample1.best")


@env.task
async def test_popscle_freemuxlet() -> None:
    # upstream case: popscle/freemuxlet "demultiplexing"
    result = await popscle_freemuxlet(await popscle_pileup(), n_sample=2, prefix="sample1")
    await check_named(result.samples, "sample1.clust1.samples.gz")
    await check_named(result.vcf, "sample1.clust1.vcf.gz")
    await check_named(result.lmix, "sample1.lmix")
    assert result.singlet_samples is None and result.singlet_vcf is None, "clust0 outputs without --aux-files"


@env.task
async def test_popscle_freemuxlet_aux_files() -> None:
    # upstream case: popscle/freemuxlet "demultiplexing - auxiliary - files"
    result = await popscle_freemuxlet(await popscle_pileup(), n_sample=2, prefix="sample1", args="--aux-files")
    await check_named(result.samples, "sample1.clust1.samples.gz")
    await check_named(result.vcf, "sample1.clust1.vcf.gz")
    await check_named(result.lmix, "sample1.lmix")
    await check_named(result.singlet_samples, "sample1.clust0.samples.gz")
    await check_named(result.singlet_vcf, "sample1.clust0.vcf.gz")


tests = [
    test_popscle_dscpileup,
    test_popscle_demuxlet_bam,
    test_popscle_demuxlet_pileup,
    test_popscle_freemuxlet,
    test_popscle_freemuxlet_aux_files,
]
