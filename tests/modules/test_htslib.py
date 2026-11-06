"""Tests for flyte_bio.modules.htslib.

md5s from upstream's htslib/bgziptabix snapshot (modules commit 9a1cc9f) and
tabix/tabix and tabix/bgzip snapshots (commit 2fb71fb, their last working
version). nf-test hashes ``.gz`` outputs after decompressing them; indexes
(``.tbi``, ``.csi``, ``.gzi``) raw.
"""

import gzip
import hashlib

from flyte.io import File

from flyte_bio.modules.htslib import htslib_bgziptabix, tabix_bgzip, tabix_tabix
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture


@env.task
async def test_bgziptabix_vcf() -> None:
    # upstream case: htslib/bgziptabix "sarscov2 - vcf - compress - index"
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf")
    gz, tbi = await htslib_bgziptabix(vcf, prefix="test", out_ext="vcf")
    assert gz.path.endswith("/test.vcf.gz") and tbi is not None, (gz.path, tbi)
    await assert_gunzipped_md5(gz, "8e722884ffb75155212a3fc053918766", label="test.vcf.gz")
    await assert_md5(tbi, "7f005943c935f2b55ba3f9d4802aa09f", label="test.vcf.gz.tbi")


@env.task
async def test_bgzip_no_index() -> None:
    # upstream case: htslib/bgziptabix "sarscov2 - gzip - (re)compress - no index" (a gzipped input)
    fastq = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    gz, tbi = await htslib_bgziptabix(fastq, prefix="test", out_ext="fastq", make_index=False)
    assert tbi is None
    await assert_gunzipped_md5(gz, "4161df271f9bfcd25d5845a1e220dbec", label="test.fastq.gz")


# tabix/tabix and tabix/bgzip: their last working upstream version (modules commit 2fb71fb).


async def variants_md5(file: File) -> str:
    """A gzipped VCF's records' md5, as nft-vcf's ``getVariantsMD5``: the non-header lines, joined without newlines."""
    async with file.open("rb") as fh:
        text = gzip.decompress(bytes(await fh.read())).decode()
    h = hashlib.md5()
    for line in text.splitlines():
        if not line.startswith("#"):
            h.update(line.encode())
    return h.hexdigest().lstrip("0")


@env.task
async def test_tabix_bed_tbi() -> None:
    # upstream case: tabix/tabix "sarscov2_bedgz_tbi"
    bed = await fixture("genomics/sarscov2/genome/bed/test.bed.gz")
    index, extracted = await tabix_tabix(bed, prefix="tbi_bed", args="-p bed")
    assert index is not None and extracted is None
    assert index.path.endswith("/test.bed.gz.tbi"), index.path
    await assert_md5(index, "9a761d51cc81835fd1199201fdbcdd5d", label="test.bed.gz.tbi")


@env.task
async def test_tabix_gff_tbi() -> None:
    # upstream case: tabix/tabix "sarscov2_gff_tbi"
    gff = await fixture("genomics/sarscov2/genome/genome.gff3.gz")
    index, _ = await tabix_tabix(gff, prefix="tbi_gff", args="-p gff")
    assert index is not None
    await assert_md5(index, "f79a67d95a98076e04fbe0455d825926", label="genome.gff3.gz.tbi")


@env.task
async def test_tabix_vcf_tbi() -> None:
    # upstream case: tabix/tabix "sarscov2_vcf_tbi"
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz")
    index, _ = await tabix_tabix(vcf, prefix="tbi_vcf", args="-p vcf")
    assert index is not None
    await assert_md5(index, "d22e5b84e4fcd18792179f72e6da702e", label="test.vcf.gz.tbi")


@env.task
async def test_tabix_vcf_csi() -> None:
    # upstream case: tabix/tabix "sarscov2_vcf_csi"
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz")
    index, _ = await tabix_tabix(vcf, prefix="vcf_csi", args="-p vcf --csi")
    assert index is not None and index.path.endswith("/test.vcf.gz.csi"), index
    await assert_md5(index, "04b41c1efd9ab3c6b1e008a286e27d2b", label="test.vcf.gz.csi")


@env.task
async def test_tabix_vcf_regions() -> None:
    # upstream case: tabix/tabix "sarscov2_vcf_tbi_regions" (its regions BED is on a fork's branch, as upstream's)
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz")
    tbi = await fixture("genomics/sarscov2/illumina/vcf/test.vcf.gz.tbi")
    bed = await fixture(
        "https://raw.githubusercontent.com/luisas/test-datasets/refs/heads/add-bedgraph-subset-illumina"
        "/data/genomics/sarscov2/illumina/bed/test.bed"
    )
    index, extracted = await tabix_tabix(vcf, tai=tbi, regions=bed, prefix="vcf_regions", args="-h")
    assert index is None and extracted is not None
    assert extracted.path.endswith("/vcf_regions.vcf.gz"), extracted.path
    actual = await variants_md5(extracted)
    assert actual == "c57539b0a39bee90f9dc6f4c60b268ed", f"variants md5 mismatch (got {actual})"


@env.task
async def test_tabix_bgzip_vcf() -> None:
    # upstream case: tabix/bgzip "sarscov2_vcf_bgzip_compress"
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf")
    gz, gzi = await tabix_bgzip(vcf, prefix="bgzip_test")
    assert gz.path.endswith("/bgzip_test.vcf.gz") and gzi is None, (gz.path, gzi)
    await assert_gunzipped_md5(gz, "8e722884ffb75155212a3fc053918766", label="bgzip_test.vcf.gz")


@env.task
async def test_tabix_bgzip_bedgz_decompress() -> None:
    # upstream case: tabix/bgzip "homo_genome_bedgz_compress" (a .gz input is decompressed)
    bed = await fixture("genomics/homo_sapiens/genome/genome.bed.gz")
    out, gzi = await tabix_bgzip(bed, prefix="bedgz_test")
    assert out.path.endswith("/bedgz_test.bed") and gzi is None, (out.path, gzi)
    await assert_md5(out, "87a15eb9c2ff20ccd5cd8735a28708f7", label="bedgz_test.bed")


@env.task
async def test_tabix_bgzip_vcf_gzi() -> None:
    # upstream case: tabix/bgzip "sarscov2_vcf_bgzip_compress_gzi" (ext.args ' -i')
    vcf = await fixture("genomics/sarscov2/illumina/vcf/test.vcf")
    gz, gzi = await tabix_bgzip(vcf, prefix="gzi_compress_test", args=" -i")
    assert gzi is not None and gzi.path.endswith("/gzi_compress_test.vcf.gz.gzi"), gzi
    await assert_gunzipped_md5(gz, "8e722884ffb75155212a3fc053918766", label="gzi_compress_test.vcf.gz")
    await assert_md5(gzi, "26fd00d4e26141cd11561f6e7d4a2ad0", label="gzi_compress_test.vcf.gz.gzi")


tests = [
    test_bgziptabix_vcf,
    test_bgzip_no_index,
    test_tabix_bed_tbi,
    test_tabix_gff_tbi,
    test_tabix_vcf_tbi,
    test_tabix_vcf_csi,
    test_tabix_vcf_regions,
    test_tabix_bgzip_vcf,
    test_tabix_bgzip_bedgz_decompress,
    test_tabix_bgzip_vcf_gzi,
]
