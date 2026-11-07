"""Tests for flyte_bio.modules.macs2.

The cases and expected md5s are upstream's macs2/callpeak tests (modules
commit 9a1cc9f, MACS2 2.2.9.1), with the ``args`` each test sets. As
upstream's inputs carry no ``single_end``, the format defaults to BAMPE; the
BED case overrides it with ``--format BED``. The ``.xls`` records MACS2's
command line, which matches upstream's because the reads keep their names.
"""

from flyte_bio.modules.macs2 import Macs2CallpeakResult, macs2_callpeak
from tests.framework import assert_md5, env, fixture

BAM = "genomics/homo_sapiens/illumina/bam/"
PACBIO_BED = (
    "genomics/homo_sapiens/pacbio/bed/"
    "alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned_tc.bed"
)


async def check(result: Macs2CallpeakResult, peak: str, xls: str, bed: str) -> None:
    assert result.peak.name == "test_peaks.narrowPeak", result.peak.name
    assert result.bed is not None and result.bed.name == "test_summits.bed"
    assert result.gapped is None and result.bdg == [], "macs2 wrote gapped peaks or pileups"
    await assert_md5(result.peak, peak, label="narrowPeak")
    await assert_md5(result.xls, xls, label="xls")
    await assert_md5(result.bed, bed, label="summits")


@env.task
async def test_macs2_callpeak_bed() -> None:
    result = await macs2_callpeak(
        [await fixture(PACBIO_BED)], gsize="4000", prefix="test", args="--format BED --qval 10 --nomodel --extsize 200"
    )
    await check(
        result,
        "10e7d4747f8a2513e5ebb04856a51673",
        "0b7c5a46179fe9d3f61c8dbc192a3c3d",
        "28833eeb7816688f0d698f51670be946",
    )


@env.task
async def test_macs2_callpeak_bam() -> None:
    result = await macs2_callpeak(
        [await fixture(BAM + "test.paired_end.name.sorted.bam")], gsize="40000", prefix="test", args="--qval 0.1"
    )
    await check(
        result,
        "2e4da1c1704595e12aaf99cc715ad70c",
        "dd0cbdd9520b150b3dd5f7bede0d4a1e",
        "26f0f97b6c14dbca129e947a58067c82",
    )


@env.task
async def test_macs2_callpeak_control_bam() -> None:
    result = await macs2_callpeak(
        [await fixture(BAM + "test.paired_end.name.sorted.bam")],
        gsize="40000",
        control=[await fixture(BAM + "test2.paired_end.name.sorted.bam")],
        prefix="test",
        args="--qval 0.1",
    )
    await check(
        result,
        "653e1108cc57ca07d0f60fc0f4fb8ba3",
        "ba5c031a290fc98828d7a3b9320863ac",
        "4f3c7c53a1d730d90d1b3dd9d3197af4",
    )


tests = [test_macs2_callpeak_bed, test_macs2_callpeak_bam, test_macs2_callpeak_control_bam]
