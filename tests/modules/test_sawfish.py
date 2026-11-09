"""Tests for flyte_bio.modules.sawfish.

The cases and expected values are upstream's sawfish/discover and
sawfish/jointcall tests (modules commit 9a1cc9f, sawfish 2.3.1), with their
configs' ``--cov-regex "."`` (unquoted here, as ``args`` isn't parsed by a
shell), and for the MAF case ``--maf-sample-name NA03697B2`` with the
sample's small variants. Discover runs on a downsampled PacBio BAM against
``genome3.fasta`` with prefix ``test``, as upstream (its contig BAM's
header records the output directory's name); its candidate BCF and index
are checked by name and its other outputs by md5, as upstream. For joint
calling, discover runs on two BAMs named ``sample1`` and ``sample2``, as
upstream's setup names them.

Joint calling runs on both samples, on ``sample1`` alone, on both read from
a sample CSV, and on both discovered with allele frequencies. As upstream,
the VCF's lines 3 to 40 are checked (here by their md5, computed from
upstream's snapshot, after the command line and sample columns are checked
directly; the first two lines hold the file date), with the index's name
and, for the MAF case, each sample's ``maf.bw`` (sawfish writes them under
``samples/sample0001_NA03697B2/`` and ``samples/sample0002_NA037562/``). The VCF header records the
command line, so the inputs keep upstream's names and the 4 CPUs its tests
ran with.
"""

import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.sawfish import SawfishDiscoverResult, SawfishJointcallResult, sawfish_discover, sawfish_jointcall
from tests.framework import assert_md5, env, fixture

BAMS = [
    "genomics/homo_sapiens/pacbio/bam/NA03697B2_downsampled.pbmm2.repeats.bam",
    "genomics/homo_sapiens/pacbio/bam/NA037562_downsampled.pbmm2.repeats.bam",
]
MAF = "genomics/homo_sapiens/pacbio/vcf/NA03697B2_new.pbmm2.repeats.vcf.gz"
FASTA = "genomics/homo_sapiens/genome/genome3.fasta"
ARGS = "--cov-regex ."
MAF_ARGS = "--cov-regex . --maf-sample-name NA03697B2"
DISCOVER_MD5 = {
    "depth.mpack": "68ca62a8de3eca85547d87735de971ac",
    "contig.alignment.bam": "a08d425c1958dcbacd9d42f9778994ef",
    "contig.alignment.bam.csi": "86e71f35df21d965390a10ce8570c125",
    "debug.breakpoint_clusters.bed": "8ed474a37f2e086c196781654a1989ad",
    "copynum.mpack": "9e16451305ea34a898fe5568ae58f852",
    "copynum.bedgraph": "2ed50e8db634ccd1a0ebc89e133b674b",
    "genome.gclevels.mpack": "f015bc227387c34d811b12bb4a81fe51",
    "sample.gcbias.mpack": "5bf4c19ef75bf38d03bc2cce769b8d50",
}
# The md5 of the VCF's lines 3 to 40 joined by newlines, from upstream's snapshot of each case.
TWO_SAMPLES_LINES_MD5 = "94734dd45b181713744c30c8db0c6db4"
ONE_SAMPLE_LINES_MD5 = "31df30f7119a398da57ab8b8297fba15"
CSV_LINES_MD5 = "0aa7e6fe8433e9dacc2c39eea9f35b85"
CMDLINE = '##sawfish_cmdline="sawfish joint-call --threads 4 --ref genome3.fasta {samples} --output-dir cohort"'


async def discover(sample: int, maf: bool = False, prefix: str = "") -> SawfishDiscoverResult:
    return await sawfish_discover(
        await fixture(BAMS[sample]),
        await fixture(BAMS[sample] + ".bai"),
        await fixture(FASTA),
        maf_vcf=await fixture(MAF) if maf else None,
        prefix=prefix or f"sample{sample + 1}",
        args=MAF_ARGS if maf else ARGS,
    )


async def jointcall(samples: int, maf: bool = False, sample_csv: File | None = None) -> SawfishJointcallResult:
    discovered = [await discover(i, maf=maf) for i in range(samples)]
    return await sawfish_jointcall(
        discovered,
        await fixture(FASTA),
        [await fixture(b) for b in BAMS[:samples]],
        [await fixture(b + ".bai") for b in BAMS[:samples]],
        sample_csv=sample_csv,
        prefix="cohort",
    )


async def check_vcf(result: SawfishJointcallResult, samples: str, columns: list[str], lines_md5: str) -> None:
    assert result.vcf.name == "cohort_genotyped.sv.vcf.gz", result.vcf.name
    assert result.tbi.name == "cohort_genotyped.sv.vcf.gz.tbi", result.tbi.name
    async with result.vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    cmdline = CMDLINE.format(samples=samples)
    assert lines[4] == cmdline, f"command line {lines[4]!r}, expected {cmdline!r}"
    chrom = [line for line in lines if line.startswith("#CHROM")]
    assert chrom and chrom[0].split("\t")[9:] == columns, chrom
    actual = hashlib.md5("\n".join(lines[2:40]).encode()).hexdigest()
    assert actual == lines_md5, f"VCF lines 3-40 md5 {actual}, expected {lines_md5}"


@env.task
async def test_sawfish_discover_bam() -> None:
    result = await discover(0, prefix="test")
    assert result.candidate_sv_bcf.name == "candidate.sv.bcf", result.candidate_sv_bcf.name
    assert result.candidate_sv_bcf_csi.name == "candidate.sv.bcf.csi", result.candidate_sv_bcf_csi.name
    for name, md5 in DISCOVER_MD5.items():
        file = await result.discover_dir.get_file(name)
        assert file is not None, f"sawfish discover wrote no {name}"
        await assert_md5(file, md5, label=name)


@env.task
async def test_sawfish_jointcall_two_samples() -> None:
    result = await jointcall(2)
    await check_vcf(result, "--sample sample1 --sample sample2", ["NA03697B2", "NA037562"], TWO_SAMPLES_LINES_MD5)


@env.task
async def test_sawfish_jointcall_one_sample() -> None:
    result = await jointcall(1)
    await check_vcf(result, "--sample sample1", ["NA03697B2"], ONE_SAMPLE_LINES_MD5)


@env.task
async def test_sawfish_jointcall_sample_csv() -> None:
    path = Path(tempfile.mkdtemp()) / "samples.csv"
    path.write_text("sample1\nsample2\n")
    result = await jointcall(2, sample_csv=await File.from_local(str(path)))
    await check_vcf(result, "--sample-csv samples.csv", ["NA03697B2", "NA037562"], CSV_LINES_MD5)


@env.task
async def test_sawfish_jointcall_maf() -> None:
    result = await jointcall(2, maf=True)
    await check_vcf(result, "--sample sample1 --sample sample2", ["NA03697B2", "NA037562"], TWO_SAMPLES_LINES_MD5)
    # sawfish names each sample's folder sample<NNNN>_<sample name>.
    maf_bw = sorted([f async for f in result.results.walk() if f.name == "maf.bw"], key=lambda f: f.path)
    folders = [f.path.split("/")[-2] for f in maf_bw]
    assert folders == ["sample0001_NA03697B2", "sample0002_NA037562"], folders
    for file, folder in zip(maf_bw, folders, strict=True):
        await assert_md5(file, "594ad937321d34420fde12ce06719efd", label=f"{folder}/maf.bw")


tests = [
    test_sawfish_discover_bam,
    test_sawfish_jointcall_two_samples,
    test_sawfish_jointcall_one_sample,
    test_sawfish_jointcall_sample_csv,
    test_sawfish_jointcall_maf,
]
