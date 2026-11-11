"""Tests for flyte_bio.modules.deepvariant.

The cases and expected md5s are upstream's deepvariant/rundeepvariant tests
(modules commit 9a1cc9f, DeepVariant 1.10.0), all with the WGS model and the
prefix ``test_out``: the homo_sapiens test BAM within ``chr22:0-40001``
against the plain and the bgzipped reference, the test CRAM within
``genome.bed``, and the CRAM again with chr22 haploid and
``genome.blacklist_intervals.bed`` as its pseudo-autosomal regions.

nf-test hashes the ``.vcf.gz`` outputs decompressed and the ``.tbi`` indexes
as they are. The VCFs, gVCFs and VCF indexes match upstream's md5s; the gVCF
indexes don't (see :func:`check`), so they get a content check.

Upstream quotes its region (``--regions="chr22:0-40001"``), which its shell
strips; ``args`` here isn't parsed by a shell, so it has no quotes.

The split steps' cases are upstream's deepvariant/makeexamples,
deepvariant/callvariants, deepvariant/postprocessvariants and
deepvariant/vcfstatsreport tests (same commit), with their configs' args
(the WGS model, or the PacBio one for the phasing case, and the small
model). make_examples' and call_variants' tests pin 2 CPUs, as upstream's
configs do, since the shard count is in every file name; the cache is
salted to keep those runs apart. The examples aren't byte-stable, so (as
upstream) only their names are checked. postprocess_variants runs on the
chain make_examples → call_variants → postprocess_variants with upstream's
4 CPUs; its VCF, gVCF and VCF index match upstream's md5s, and its gVCF
index gets a content check, as above. The visual report matches upstream's
md5 for VCF, bgzipped VCF, BCF and compressed BCF input.
"""

import gzip

import flyte
from flyte.io import File

from flyte_bio.modules.deepvariant import (
    DeepVariantExamples,
    DeepVariantResult,
    deepvariant_callvariants,
    deepvariant_callvariants_cmd,
    deepvariant_examples,
    deepvariant_makeexamples,
    deepvariant_makeexamples_cmd,
    deepvariant_postprocessvariants,
    deepvariant_rundeepvariant,
    deepvariant_vcfstatsreport,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
BAM = "genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam"
CRAM = "genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram"

DIPLOID = {
    "vcf": "d6ec6971b53fc10bf9ed23ddf6f544c1",
    "vcf_index": "efb861c33cf51739a4ae13a2e81126b6",
    "gvcf": "71bf2757abd13abed78479887f9de8f6",
    "gvcf_index": "90732f1ea2188c12fc4dab62da039f63",  # not reproduced, see check()
}
HAPLOID_CHR22 = {
    "vcf": "44664d7a287819aa3c18d5bd1aafed0b",
    "vcf_index": "7f57f737784a7102360ae06a636eac14",
    "gvcf": "61246dde2499b2301429eff9e96714dd",
    "gvcf_index": "3cb9cc54bbcd6bf44d2b835aeeaa5f9e",  # not reproduced, see check()
}


async def check(result: DeepVariantResult, expected: dict[str, str]) -> None:
    """Check a VCF and gVCF named ``test_out.*`` against ``expected``."""
    await check_named(result, expected, "test_out")


async def check_named(result: DeepVariantResult, expected: dict[str, str], prefix: str) -> None:
    files: dict[str, File] = {
        "vcf": result.vcf,
        "vcf_index": result.vcf_index,
        "gvcf": result.gvcf,
        "gvcf_index": result.gvcf_index,
    }
    names = {key: file.name for key, file in files.items()}
    assert names == {
        "vcf": f"{prefix}.vcf.gz",
        "vcf_index": f"{prefix}.vcf.gz.tbi",
        "gvcf": f"{prefix}.g.vcf.gz",
        "gvcf_index": f"{prefix}.g.vcf.gz.tbi",
    }, names
    for key in ("vcf", "vcf_index", "gvcf"):
        assert_file_md5 = assert_gunzipped_md5 if key in ("vcf", "gvcf") else assert_md5
        await assert_file_md5(files[key], expected[key], label=files[key].name)
    # The gVCF's index records the gVCF's compressed block offsets, and our
    # gVCF compresses to different bytes than upstream's (its content matches):
    # check it's a tabix index instead of its md5.
    async with result.gvcf_index.open("rb") as fh:
        index = gzip.decompress(bytes(await fh.read()))
    assert index[:4] == b"TBI\x01", f"{result.gvcf_index.name} isn't a tabix index"
    assert result.report is None, result.report


@env.task
async def test_deepvariant_rundeepvariant_bam() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(BAM),
        await fixture(BAM + ".bai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        prefix="test_out",
        args="--regions=chr22:0-40001 --model_type=WGS",
    )
    await check(result, DIPLOID)


@env.task
async def test_deepvariant_rundeepvariant_cram_intervals() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(CRAM),
        await fixture(CRAM + ".crai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        intervals=await fixture(GENOME + "genome.bed"),
        prefix="test_out",
        args="--model_type=WGS",
    )
    await check(result, DIPLOID)


@env.task
async def test_deepvariant_rundeepvariant_cram_intervals_par_bed() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(CRAM),
        await fixture(CRAM + ".crai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        intervals=await fixture(GENOME + "genome.bed"),
        par_bed=await fixture(GENOME + "genome.blacklist_intervals.bed"),
        prefix="test_out",
        args="--model_type=WGS --haploid_contigs chr22",
    )
    await check(result, HAPLOID_CHR22)


@env.task
async def test_deepvariant_rundeepvariant_bam_fasta_gz() -> None:
    result = await deepvariant_rundeepvariant(
        await fixture(BAM),
        await fixture(BAM + ".bai"),
        await fixture(GENOME + "genome.fasta.gz"),
        await fixture(GENOME + "genome.fasta.gz.fai"),
        gzi=await fixture(GENOME + "genome.fasta.gz.gzi"),
        prefix="test_out",
        args="--regions=chr22:0-40001 --model_type=WGS",
    )
    await check(result, DIPLOID)


# The split steps' cases: upstream's test configs' args, unquoted.
WGS_EXAMPLES_ARGS = (
    "--checkpoint /opt/models/wgs --call_small_model_examples --small_model_indel_gq_threshold 30"
    " --small_model_snp_gq_threshold 25 --small_model_vaf_context_window_size 51"
    " --trained_small_model_path /opt/smallmodels/wgs/model.keras"
)
PACBIO_EXAMPLES_ARGS = (
    "--checkpoint /opt/models/pacbio --call_small_model_examples --small_model_indel_gq_threshold 30"
    " --small_model_snp_gq_threshold 25 --small_model_vaf_context_window_size 51"
    " --trained_small_model_path /opt/smallmodels/pacbio/model.keras"
)
EXAMPLES_NAMES = [
    "test.examples.tfrecord-00000-of-00002.gz",
    "test.examples.tfrecord-00000-of-00002.gz.example_info.json",
    "test.examples.tfrecord-00001-of-00002.gz",
    "test.examples.tfrecord-00001-of-00002.gz.example_info.json",
]
GVCF_NAMES = ["test.gvcf.tfrecord-00000-of-00002.gz", "test.gvcf.tfrecord-00001-of-00002.gz"]
POSTPROCESS = {
    "vcf": "d6ec6971b53fc10bf9ed23ddf6f544c1",
    "vcf_index": "efb861c33cf51739a4ae13a2e81126b6",
    "gvcf": "71bf2757abd13abed78479887f9de8f6",
    "gvcf_index": "7a0dae7f0b14f5f62b00d5edf0f5d4f7",  # not reproduced, see check()
}
# Upstream's tests' CPUs for make_examples and call_variants; the salt keeps these runs' cache entries apart
# from default-resource ones.
TWO_CPUS = {
    "resources": flyte.Resources(cpu=2, memory="16Gi", disk="20Gi"),
    "cache": flyte.Cache(behavior="auto", salt="cpu=2"),
}


async def deepvariant_makeexamples_two_cpus(
    reads: File,
    index: File,
    fasta: File,
    fai: File,
    intervals: File | None = None,
    gzi: File | None = None,
    with_phasing: bool = False,
    args: str = WGS_EXAMPLES_ARGS,
) -> DeepVariantExamples:
    task = deepvariant_makeexamples_cmd.override(**TWO_CPUS)
    results = await task(
        reads=reads,
        index=index,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        par_bed=None,
        with_phasing=with_phasing,
        prefix="test",
        args=args,
    )
    return await deepvariant_examples(results, "test")


async def bam_examples_two_cpus(with_phasing: bool = False, args: str = WGS_EXAMPLES_ARGS) -> DeepVariantExamples:
    return await deepvariant_makeexamples_two_cpus(
        await fixture(BAM),
        await fixture(BAM + ".bai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        with_phasing=with_phasing,
        args=args,
    )


def check_examples(examples: DeepVariantExamples, phasing: bool = False) -> None:
    names = [f.name for f in examples.examples]
    assert names == EXAMPLES_NAMES, names
    names = [f.name for f in examples.gvcf]
    assert names == GVCF_NAMES, names
    names = [f.name for f in examples.read_phase_inputs]
    expected = ["test-read-phasing_debug-00000-of-00002.tsv", "test-read-phasing_debug-00001-of-00002.tsv"]
    assert names == (expected if phasing else []), names


@env.task
async def test_deepvariant_makeexamples_bam() -> None:
    check_examples(await bam_examples_two_cpus())


@env.task
async def test_deepvariant_makeexamples_cram_intervals() -> None:
    examples = await deepvariant_makeexamples_two_cpus(
        await fixture(CRAM),
        await fixture(CRAM + ".crai"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.fasta.fai"),
        intervals=await fixture(GENOME + "genome.bed"),
    )
    check_examples(examples)


@env.task
async def test_deepvariant_makeexamples_bam_fasta_gz() -> None:
    examples = await deepvariant_makeexamples_two_cpus(
        await fixture(BAM),
        await fixture(BAM + ".bai"),
        await fixture(GENOME + "genome.fasta.gz"),
        await fixture(GENOME + "genome.fasta.gz.fai"),
        gzi=await fixture(GENOME + "genome.fasta.gz.gzi"),
    )
    check_examples(examples)


@env.task
async def test_deepvariant_makeexamples_bam_phasing() -> None:
    check_examples(await bam_examples_two_cpus(with_phasing=True, args=PACBIO_EXAMPLES_ARGS), phasing=True)


@env.task
async def test_deepvariant_callvariants_wgs() -> None:
    examples = await bam_examples_two_cpus()
    task = deepvariant_callvariants_cmd.override(**TWO_CPUS)
    results = await task(examples=examples.examples, prefix="test", args="--checkpoint /opt/models/wgs")
    names = sorted([f.name async for f in results.walk()])
    assert names == ["test.call-00000-of-00001.tfrecord.gz"], names


@env.task
async def test_deepvariant_postprocessvariants_wgs() -> None:
    fasta = await fixture(GENOME + "genome.fasta")
    fai = await fixture(GENOME + "genome.fasta.fai")
    examples = await deepvariant_makeexamples(
        await fixture(BAM), await fixture(BAM + ".bai"), fasta, fai, prefix="test", args=WGS_EXAMPLES_ARGS
    )
    calls = await deepvariant_callvariants(examples.examples, prefix="test", args="--checkpoint /opt/models/wgs")
    result = await deepvariant_postprocessvariants(
        calls,
        examples.gvcf,
        fasta,
        fai,
        small_model_calls=examples.small_model_calls,
        phase_inputs=examples.read_phase_inputs,
        prefix="test",
    )
    await check_named(result, POSTPROCESS, "test")


async def check_report(path: str) -> None:
    report = await deepvariant_vcfstatsreport(await fixture(path), prefix="test")
    assert report.name == "test.visual_report.html", report.name
    await assert_md5(report, "7b6b41bc83da2e608416f2572f0e9f9d", label=report.name)


@env.task
async def test_deepvariant_vcfstatsreport_vcf() -> None:
    await check_report("genomics/sarscov2/illumina/vcf/test.vcf")


@env.task
async def test_deepvariant_vcfstatsreport_vcf_gz() -> None:
    await check_report("genomics/sarscov2/illumina/vcf/test.vcf.gz")


@env.task
async def test_deepvariant_vcfstatsreport_bcf() -> None:
    await check_report("genomics/sarscov2/illumina/vcf/test.bcf")


@env.task
async def test_deepvariant_vcfstatsreport_bcf_gz() -> None:
    await check_report("genomics/sarscov2/illumina/vcf/test.bcf.gz")


tests = [
    test_deepvariant_rundeepvariant_bam,
    test_deepvariant_rundeepvariant_cram_intervals,
    test_deepvariant_rundeepvariant_cram_intervals_par_bed,
    test_deepvariant_rundeepvariant_bam_fasta_gz,
    test_deepvariant_makeexamples_bam,
    test_deepvariant_makeexamples_cram_intervals,
    test_deepvariant_makeexamples_bam_fasta_gz,
    test_deepvariant_makeexamples_bam_phasing,
    test_deepvariant_callvariants_wgs,
    test_deepvariant_postprocessvariants_wgs,
    test_deepvariant_vcfstatsreport_vcf,
    test_deepvariant_vcfstatsreport_vcf_gz,
    test_deepvariant_vcfstatsreport_bcf,
    test_deepvariant_vcfstatsreport_bcf_gz,
]
