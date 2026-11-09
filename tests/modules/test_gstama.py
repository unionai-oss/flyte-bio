"""Tests for flyte_bio.modules.gstama.

The cases and expected md5s are upstream's gstama/collapse, gstama/merge and
gstama/polyacleanup tests (modules commit 9a1cc9f, gs-tama 1.0.3). Collapse
runs with its test config's ``-x capped -b BAM`` and prefix ``<id>_tc``.
nf-test hashes every output, the gzipped ones decompressed.
"""

import asyncio

from flyte_bio.modules.gstama import gstama_collapse, gstama_merge, gstama_polyacleanup
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

PACBIO = "genomics/homo_sapiens/pacbio/"
ALIGNED = "alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.singletons.merged.aligned"


@env.task
async def test_gstama_collapse() -> None:
    # upstream case: gstama/collapse "test-gstama-collapse"
    bam, fasta = await asyncio.gather(
        fixture(PACBIO + f"bam/{ALIGNED}.bam"), fixture("genomics/homo_sapiens/genome/genome2.fasta")
    )
    result = await gstama_collapse(bam, fasta, prefix="test_tc", args="-x capped -b BAM")
    assert result.varcov is not None and result.variants is not None, "tama_collapse.py called no variants"
    expected = {
        "bed": ("test_tc_collapsed.bed", "e5105198ed970a33ae0ecaa7bff421d9"),
        "bed_trans_reads": ("test_tc_trans_read.bed", "0ca1a32f33ef05242d897d913802554b"),
        "local_density_error": ("test_tc_local_density_error.txt", "b917ac1f14eccd590b6881a686f324d5"),
        "polya": ("test_tc_polya.txt", "628ea62b918fc4f31e109f724d714a66"),
        "read": ("test_tc_read.txt", "d2685d7f24cd1611e0770a5ce25422fe"),
        "strand_check": ("test_tc_strand_check.txt", "42cc52b2660b1e0b84e1c9ab37a965ec"),
        "trans_report": ("test_tc_trans_report.txt", "33a86c15ca2acce36b2a5962f4c1adc4"),
        "varcov": ("test_tc_varcov.txt", "587fd899ff658eb66b1770a35283bfcb"),
        "variants": ("test_tc_variants.txt", "5b1165e9f33faba4f7207013fc27257e"),
    }
    for key, (name, md5) in expected.items():
        await assert_md5(getattr(result, key), md5, label=name)


@env.task
async def test_gstama_merge() -> None:
    # upstream case: gstama/merge "test-gstama-merge"
    beds = list(
        await asyncio.gather(
            fixture(PACBIO + f"bed/{ALIGNED}_tc.bed"),
            fixture(PACBIO + f"bed/{ALIGNED}_tc.2.bed"),
        )
    )
    filelist = await fixture(PACBIO + "txt/filelist.txt")
    result = await gstama_merge(beds, filelist, prefix="test_all")
    assert result.bed.name == "test_all.bed", result.bed.name
    await assert_md5(result.bed, "60ec34e1ff9655d4ce2e83d3f4bbf448", label="test_all.bed")
    await assert_md5(result.gene_report, "7029fd183dfd905a233403cfbe44722a", label="test_all_gene_report.txt")
    await assert_md5(result.merge, "4279e59ed5739ce4f2f811568962893f", label="test_all_merge.txt")
    await assert_md5(result.trans_report, "97d8346d9eb9da140941656c3a3325cd", label="test_all_trans_report.txt")


@env.task
async def test_gstama_polyacleanup() -> None:
    # upstream case: gstama/polyacleanup "test-gstama-polyacleanup"
    fasta = await fixture("genomics/homo_sapiens/genome/transcriptome.fasta")
    result = await gstama_polyacleanup(fasta, prefix="test")
    assert result.fasta.name == "test.fa.gz", result.fasta.name
    await assert_gunzipped_md5(result.fasta, "f0a2c3ca8f19d2197c6b5b273093ebf0", label="test.fa.gz")
    await assert_gunzipped_md5(result.report, "83ba9d11a59ff516f45dd45c78613206", label="test_polya_flnc_report.txt.gz")
    await assert_gunzipped_md5(result.tails, "20a02db3d29cc45b39880e1ad5ee243b", label="test_tails.fa.gz")


tests = [test_gstama_collapse, test_gstama_merge, test_gstama_polyacleanup]
