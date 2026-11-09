"""Tests for flyte_bio.modules.vcf2maf.

The cases are upstream's vcf2maf tests (modules commit 9a1cc9f, vcf2maf
1.6.22, VEP 113), with their ``ext.args``:

- without VEP (``--inhibit-vep``): as upstream, the MAF's lines 6 to 132 (its
  first 127 variant rows) are compared; their md5 over ``"\\n".join(lines)``
  is computed from upstream's snapshotted lines;
- with VEP and upstream's tiny test cache (``vep_cache_113.tar.gz``, which
  annotates nothing): the MAF is just its two header lines, checked by md5.

Upstream also snapshots the tools' versions (vcf2maf 1.6.22, a constant, and
VEP's ``113.0`` from ``vep --help``); the task doesn't report versions.
"""

import hashlib

from flyte.io import File

from flyte_bio.modules.vcf2maf import vcf2maf
from tests.framework import assert_md5, env, fixture, fixture_dir

VCF = "genomics/homo_sapiens/illumina/gvcf/test.genome.vcf"
FASTA = "genomics/homo_sapiens/genome/genome.fasta"
# md5 of "\n".join(lines[5:132]) of upstream's "Run without VEP" snapshot.
NO_VEP_LINES_MD5 = "032b3163ac668f47a2e2c08d75e462fb"


async def lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


@env.task
async def test_vcf2maf_without_vep() -> None:
    # upstream case: "Run without VEP"
    args = "--ncbi-build GRCh38 --species homo_sapiens --inhibit-vep"
    result = await vcf2maf(await fixture(VCF), await fixture(FASTA), prefix="test", args=args)
    rows = (await lines(result.maf))[5:132]
    assert len(rows) == 127, f"test.maf has {len(rows) + 5} lines, expected at least 132"
    actual = hashlib.md5("\n".join(rows).encode()).hexdigest()
    assert actual == NO_VEP_LINES_MD5, f"test.maf lines 6-132: md5 {actual}, expected {NO_VEP_LINES_MD5}"


@env.task
async def test_vcf2maf_with_vep() -> None:
    # upstream case: "Run with VEP"
    cache = await fixture_dir("genomics/homo_sapiens/genome/vep_cache_113.tar.gz", "vep")
    args = "--ncbi-build GRCh38 --species homo_sapiens"
    result = await vcf2maf(await fixture(VCF), await fixture(FASTA), vep_cache=cache, prefix="test", args=args)
    assert len(await lines(result.maf)) == 2, "test.maf has variants; expected only its two header lines"
    await assert_md5(result.maf, "8491ed276136b31ed9e94478eb9ddf2c", label="test.maf")


tests = [test_vcf2maf_without_vep, test_vcf2maf_with_vep]
