"""Tests for flyte_bio.modules.artic.

The cases are upstream's artic/guppyplex and artic/minion tests (modules
commit 9a1cc9f, artic 1.6.2). guppyplex's output is checked by its
decompressed md5; minion's as upstream: the names of everything it wrote,
the BAMs' names, and the md5s of the BAM indexes, the consensus and the
passing VCF's index (no JSON). minion uses the image's bundled ``hifi``
Clair3 model and the ARTIC V3 SARS-CoV-2 scheme, as upstream.
"""

import tempfile
from pathlib import Path

from flyte.io import Dir

from flyte_bio.modules.artic import artic_guppyplex, artic_minion
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

NANOPORE_FASTQ = "genomics/sarscov2/nanopore/fastq/test.fastq.gz"
SCHEME = "https://github.com/artic-network/primer-schemes/raw/master/nCoV-2019/V3/"
BAI_MD5 = "0d76977b2e36046cc176112776c5fa4e"
MINION_FILES = [
    "test.1.trimmed.rg.sorted.bam",
    "test.1.trimmed.rg.sorted.bam.bai",
    "test.1.vcf",
    "test.2.trimmed.rg.sorted.bam",
    "test.2.trimmed.rg.sorted.bam.bai",
    "test.2.vcf",
    "test.alignreport.csv",
    "test.amplicon_depths.tsv",
    "test.consensus.fasta",
    "test.coverage_mask.txt",
    "test.coverage_mask.txt.1.depths",
    "test.coverage_mask.txt.2.depths",
    "test.fail.vcf",
    "test.merged.vcf",
    "test.merged.vcf.gz",
    "test.merged.vcf.gz.tbi",
    "test.minion.log.txt",
    "test.normalised.vcf.gz",
    "test.normalised.vcf.gz.tbi",
    "test.pass.vcf",
    "test.pass.vcf.gz",
    "test.pass.vcf.gz.tbi",
    "test.preconsensus.fasta",
    "test.preconsensus.fasta.fai",
    "test.primers.vcf",
    "test.primersitereport.txt",
    "test.primertrimmed.rg.sorted.bam",
    "test.primertrimmed.rg.sorted.bam.bai",
    "test.sorted.bam",
    "test.sorted.bam.bai",
    "test.trimmed.rg.sorted.bam",
    "test.trimmed.rg.sorted.bam.bai",
]


@env.task
async def test_artic_guppyplex() -> None:
    # upstream case: artic/guppyplex "sarscov2 - fastq"
    staged = Path(tempfile.mkdtemp()) / "fastq"
    staged.mkdir()
    await (await fixture(NANOPORE_FASTQ)).download(str(staged / "test.fastq.gz"))
    result = await artic_guppyplex(await Dir.from_local(str(staged)), prefix="test")
    assert result.name == "test.fastq.gz", result.name
    await assert_gunzipped_md5(result, "7567d853ada6ac142332619d0b541d76", label="test.fastq.gz")


@env.task
async def test_artic_minion() -> None:
    # upstream case: artic/minion "sarscov2 - fastq"
    result = await artic_minion(
        await fixture(NANOPORE_FASTQ),
        await fixture(SCHEME + "nCoV-2019.reference.fasta"),
        await fixture(SCHEME + "nCoV-2019.primer.bed"),
        model="hifi",
        prefix="test",
    )
    names = sorted(f.name for f in await result.results.list_files())
    assert names == MINION_FILES, f"artic minion wrote {names}"
    assert (result.bam.name, result.bam_trimmed.name, result.bam_primertrimmed.name) == (
        "test.sorted.bam",
        "test.trimmed.rg.sorted.bam",
        "test.primertrimmed.rg.sorted.bam",
    )
    for bai in (result.bai, result.bai_trimmed, result.bai_primertrimmed):
        await assert_md5(bai, BAI_MD5, label=bai.name)
    await assert_md5(result.fasta, "ab9bbc82f542c9eed46808c640aa7a1e", label="test.consensus.fasta")
    assert result.vcf.name == "test.pass.vcf.gz", result.vcf.name
    await assert_md5(result.tbi, "4cb176febbc8c26d717a6c6e67b9c905", label="test.pass.vcf.gz.tbi")
    assert result.json == [], [f.name for f in result.json]


tests = [test_artic_guppyplex, test_artic_minion]
