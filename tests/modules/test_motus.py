"""Tests for flyte_bio.modules.motus.

The cases are upstream's motus/preplong, motus/profile and motus/merge tests
(modules commit 9a1cc9f, mOTUs 3.1.0). As upstream's setups, the database is
downloaded with mOTUs' ``downloadDB.py`` (from the mOTUs repository, as
upstream; about 3 GB, cached after the first run). The prep_long output is
checked by md5 (decompressed, as nf-test's); as upstream, the profiles, logs
and merged table are checked by name (their headers record paths).
"""

from flyte.io import File

from flyte_bio.modules.motus import motus_downloaddb, motus_merge, motus_preplong, motus_profile
from tests.framework import assert_gunzipped_md5, env, fixture

DOWNLOAD_SCRIPT = "https://raw.githubusercontent.com/motu-tool/mOTUs/master/motus/downloadDB.py"
FASTQ = "genomics/sarscov2/illumina/fastq/"


async def motus_test_db() -> File:
    return await motus_downloaddb(await fixture(DOWNLOAD_SCRIPT))


@env.task
async def test_motus_preplong() -> None:
    # upstream case: motus/preplong "sarscov2 - fasta"
    reads = await fixture("genomics/sarscov2/nanopore/fastq/test_2.fastq.gz")
    converted = await motus_preplong(reads, await motus_test_db(), prefix="test")
    assert converted.name == "test.gz", converted.name
    await assert_gunzipped_md5(converted, "37be4a7321bcb2fdae60affb21a1737c", label="test.gz")


@env.task
async def test_motus_profile_fastq() -> None:
    # upstream case: motus/profile "sarscov2 - illumina paired end [fastq]"
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    result = await motus_profile(reads, await motus_test_db(), prefix="test")
    assert (result.out.name, result.log.name) == ("test.out", "test.log"), (result.out.name, result.log.name)
    assert result.bam is None and result.mgc is None


@env.task
async def test_motus_profile_bam() -> None:
    # upstream case: motus/profile "sarscov2 - illumina paired end [bam]"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam")
    result = await motus_profile(bam, await motus_test_db(), prefix="test")
    assert (result.out.name, result.log.name) == ("test.out", "test.log"), (result.out.name, result.log.name)
    assert result.bam is None and result.mgc is None


@env.task
async def test_motus_merge() -> None:
    # upstream case: motus/merge "sarscov2 - illumina paired end [fastq]"
    db = await motus_test_db()
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    profile = await motus_profile(reads, db, prefix="test")
    merged = await motus_merge([profile.out], db, prefix="test")
    assert merged.name == "test.txt", merged.name


tests = [test_motus_preplong, test_motus_profile_fastq, test_motus_profile_bam, test_motus_merge]
