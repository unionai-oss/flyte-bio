"""Tests for flyte_bio.modules.hiphase.

The cases and expected md5s are upstream's hiphase tests (modules commit
9a1cc9f, HiPhase 1.6.0): a downsampled PacBio BAM and the sample's small
variants (with their ``.csi``) against ``genome3.fasta`` (without a
``.fai``), phasing the VCF only, then also writing the haplotagged BAM and
the summary, blocks, stats and haplotag files as TSV, then as CSV. nf-test
hashes the phased ``.vcf.gz`` decompressed and the indexes and BAM as they
are; as upstream, only the report files' names are checked. HiPhase records
its command line in the VCF and BAM headers, so the module runs it with
upstream's 4 CPUs and the inputs under their own names.
"""

from flyte_bio.modules.hiphase import HiphaseResult, hiphase
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

PACBIO = "genomics/homo_sapiens/pacbio/"
BAM = PACBIO + "bam/NA03697B2_downsampled.pbmm2.repeats.bam"
VCF = PACBIO + "vcf/NA03697B2_new.pbmm2.repeats.vcf.gz"
FASTA = "genomics/homo_sapiens/genome/genome3.fasta"
HAPLOTAGGED = "NA03697B2_downsampled.pbmm2.repeats_haplotagged.bam"
VCF_INDEX_MD5 = "f4ac1b73f3279823944a5c821a709e45"
BAM_INDEX_MD5 = "c2545d124e5398d397670758f5e46fc4"


async def run(all_outputs: bool, file_format: str = "tsv") -> HiphaseResult:
    return await hiphase(
        [await fixture(BAM)],
        [await fixture(BAM + ".bai")],
        await fixture(FASTA),
        snv=await fixture(VCF),
        snv_index=await fixture(VCF + ".csi"),
        output_bam=all_outputs,
        summary_file=all_outputs,
        blocks_file=all_outputs,
        stats_file=all_outputs,
        haplotag_file=all_outputs,
        file_format=file_format,
        prefix="test",
    )


async def check_vcf(result: HiphaseResult, vcf_md5: str) -> None:
    assert result.snv_vcf is not None and result.snv_vcf.name == "test_snv_phased.vcf.gz", result.snv_vcf
    assert result.snv_vcf_index is not None, "no index of the phased VCF"
    assert result.snv_vcf_index.name == "test_snv_phased.vcf.gz.tbi", result.snv_vcf_index.name
    assert result.sv_vcf is None and result.sv_vcf_index is None
    await assert_gunzipped_md5(result.snv_vcf, vcf_md5, label=result.snv_vcf.name)
    await assert_md5(result.snv_vcf_index, VCF_INDEX_MD5, label=result.snv_vcf_index.name)


async def check_all_outputs(result: HiphaseResult, file_format: str, vcf_md5: str, bam_md5: str) -> None:
    await check_vcf(result, vcf_md5)
    reports = [result.summary, result.blocks, result.stats, result.haplotag]
    names = [r.name if r else None for r in reports]
    assert names == [f"test.{kind}.{file_format}" for kind in ("summary", "blocks", "stats", "haplotag")], names
    assert [b.name for b in result.bams] == [HAPLOTAGGED], result.bams
    assert [b.name for b in result.bam_indexes] == [HAPLOTAGGED + ".bai"], result.bam_indexes
    await assert_md5(result.bams[0], bam_md5, label=HAPLOTAGGED)
    await assert_md5(result.bam_indexes[0], BAM_INDEX_MD5, label=HAPLOTAGGED + ".bai")


@env.task
async def test_hiphase_snv_only() -> None:
    result = await run(all_outputs=False)
    await check_vcf(result, "dc898f8a5526369e84ea34bbf820ee78")
    assert result.bams == [] and result.summary is None and result.haplotag is None


@env.task
async def test_hiphase_all_outputs_tsv() -> None:
    result = await run(all_outputs=True, file_format="tsv")
    await check_all_outputs(result, "tsv", "abec5b169e3ecf1b345dbbb4de93dd07", "c18a4030683f2ad7d4889922142bbb03")


@env.task
async def test_hiphase_all_outputs_csv() -> None:
    result = await run(all_outputs=True, file_format="csv")
    await check_all_outputs(result, "csv", "b7c30e9a2b2483847e7d4dcf87a5084c", "b7ebb05ee7914006481278b54f9e8202")


tests = [
    test_hiphase_snv_only,
    test_hiphase_all_outputs_tsv,
    test_hiphase_all_outputs_csv,
]
