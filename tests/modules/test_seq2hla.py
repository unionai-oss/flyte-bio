"""Tests for flyte_bio.modules.seq2hla (md5s from upstream's seq2hla module snapshot)."""

import asyncio

from flyte_bio.modules.seq2hla import seq2hla
from tests.framework import assert_md5, env, fixture

FASTQ = "genomics/homo_sapiens/illumina/fastq/"
BOWTIELOG_MD5 = "6e932ca10d185c4f226a20a4948ab175"


@env.task
async def test_seq2hla_fastq() -> None:
    # upstream case: seq2hla "homo_sapiens - fastq". Upstream snapshots only the expression files' names
    # (on this small test sample they are empty).
    r1, r2 = await asyncio.gather(fixture(FASTQ + "test_rnaseq_1.fastq.gz"), fixture(FASTQ + "test_rnaseq_2.fastq.gz"))
    result = await seq2hla(r1, r2, prefix="test")
    expected = {
        "ClassI-class": (result.class1, "6b501035e8d1d51ca92fe37fb1c6dda7", "cc22f9df7ccdb7a50c1a16b3621422f4"),
        "ClassI-nonclass": (
            result.class1_nonclass,
            "b97d79e23ef66f946478ea913dc0e641",
            "63958dcdcd2d912b44ee9bb3ca651b10",
        ),
        "ClassII": (result.class2, "2a8996f46f7749486b2bbb3bbf15032f", "d48e48f6fb7bd4309e1071b78acaaffa"),
    }
    for hla_class, (files, md5_2d, md5_4d) in expected.items():
        stem = f"test-{hla_class}"
        assert files.expression.name == f"{stem}.expression", files.expression.name
        await assert_md5(files.bowtielog, BOWTIELOG_MD5, label=f"{stem}.bowtielog")
        await assert_md5(files.genotype_2d, md5_2d, label=f"{stem}.HLAgenotype2digits")
        await assert_md5(files.genotype_4d, md5_4d, label=f"{stem}.HLAgenotype4digits")
    assert result.ambiguity is None, "upstream writes no .ambiguity for this sample"


tests = [test_seq2hla_fastq]
