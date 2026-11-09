"""Tests for flyte_bio.modules.hifiasm.

The cases and expected md5s are upstream's hifiasm tests (modules commit
9a1cc9f, hifiasm 0.25.0). The module's 4 CPUs match upstream's test CPUs.

Upstream's three failing cases (trio with Hi-C, one k-mer dump, one Hi-C
mate) fail there because their test inputs are malformed (the module only
logs an error); here the wrapper rejects those combinations, so the tests
check that it raises.
"""

from flyte.io import File

from flyte_bio.modules.hifiasm import HifiasmResult, hifiasm
from tests.framework import assert_md5, assert_nonempty, env, fixture

HIFI = "genomics/homo_sapiens/pacbio/fastq/test_hifi.fastq.gz"
YAK = "genomics/homo_sapiens/illumina/yak/"
HIC = "genomics/homo_sapiens/illumina/fastq/"
BASE_BINS = ["test.ec.bin", "test.ovlp.reverse.bin", "test.ovlp.source.bin"]


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: HifiasmResult, expected: dict[str, tuple[str, str] | None]) -> None:
    """``expected`` maps each output to (name, md5), or None when it must be absent."""
    for label, want in expected.items():
        file = getattr(result, label)
        if want is None:
            assert file is None, f"hifiasm wrote unexpected {label}: {file.name}"
            continue
        assert file is not None, f"hifiasm wrote no {label}"
        name, md5 = want
        assert file.name == name, f"{label} is {file.name}, expected {name}"
        await assert_md5(file, md5, label=name)
    await assert_nonempty(result.log, label="log")


def check_bins(result: HifiasmResult, expected: list[str]) -> None:
    names = sorted(f.name for f in result.bin_files)
    assert names == sorted(expected), f"bin files {names}, expected {sorted(expected)}"


BP = {
    "raw_unitigs": ("test.bp.r_utg.gfa", "6a200e4bcb5d8d8fa6ed75874d26edc6"),
    "processed_unitigs": ("test.bp.p_utg.gfa", "6a200e4bcb5d8d8fa6ed75874d26edc6"),
    "primary_contigs": ("test.bp.p_ctg.gfa", "c3225425a73a8d6d46536dab72380650"),
    "hap1_contigs": ("test.bp.hap1.p_ctg.gfa", "78a86b6d5741a680e22a427ca9faa74a"),
    "hap2_contigs": ("test.bp.hap2.p_ctg.gfa", "ac2116fd2f22c67d4c304cbf9b9f7793"),
}


@env.task
async def test_hifiasm_hifi() -> None:
    result = await hifiasm([await fixture(HIFI)], prefix="test", args="-f0 --write-ec --write-paf")
    check_bins(result, BASE_BINS)
    await check(result, BP)
    assert result.corrected_reads is not None and result.read_overlaps is not None, "no --write-ec/--write-paf output"


@env.task
async def test_hifiasm_hifi_bins() -> None:
    reads = [await fixture(HIFI)]
    initial = await hifiasm(reads, prefix="test", args="-f0 --bin-only")
    result = await hifiasm(reads, bin_files=initial.bin_files, prefix="test", args="-f0")
    first = (await text(result.log)).splitlines()[0]
    assert "Reads has been loaded." in first, f"log starts {first!r}, not with the bins being loaded"
    await check(result, BP)


@env.task
async def test_hifiasm_hifi_two_files() -> None:
    reads = [await fixture(HIFI), await fixture("genomics/homo_sapiens/pacbio/fastq/alz.ccs.fastq")]
    result = await hifiasm(reads, prefix="test", args="-f0")
    check_bins(result, BASE_BINS)
    await check(
        result,
        {
            "raw_unitigs": ("test.bp.r_utg.gfa", "e6e38ac76f6b73142c3fe403f71d27b0"),
            "processed_unitigs": ("test.bp.p_utg.gfa", "e6e38ac76f6b73142c3fe403f71d27b0"),
            "primary_contigs": ("test.bp.p_ctg.gfa", "43890a1832d8f26de263e57dc5e3b8de"),
            "hap1_contigs": ("test.bp.hap1.p_ctg.gfa", "7d7ea2bed472de263f6ec3521959b0d9"),
            "hap2_contigs": ("test.bp.hap2.p_ctg.gfa", "ce096a66c9bba039c6a22ba9e9409d01"),
        },
    )


@env.task
async def test_hifiasm_hifi_trio() -> None:
    result = await hifiasm(
        [await fixture(HIFI)],
        paternal_kmer_dump=await fixture(YAK + "test.yak"),
        maternal_kmer_dump=await fixture(YAK + "test2.yak"),
        prefix="test",
        args="-f0",
    )
    check_bins(result, [*BASE_BINS, "test.hap1.phase.bin", "test.hap2.phase.bin"])
    assert result.hap2_contigs is not None and result.hap2_contigs.name == "test.dip.hap2.p_ctg.gfa"
    await check(
        result,
        {
            "raw_unitigs": ("test.dip.r_utg.gfa", "68361ac3e8babd51f3891d1637ca0fdc"),
            "processed_unitigs": ("test.dip.p_utg.gfa", "68361ac3e8babd51f3891d1637ca0fdc"),
            "hap1_contigs": ("test.dip.hap1.p_ctg.gfa", "eed5da5f3dd415dbb711edb61a09802f"),
        },
    )


@env.task
async def test_hifiasm_hifi_hic() -> None:
    result = await hifiasm(
        [await fixture(HIFI)],
        hic_read1=await fixture(HIC + "test_1.fastq.gz"),
        hic_read2=await fixture(HIC + "test_2.fastq.gz"),
        prefix="test",
        args="-f0",
    )
    check_bins(result, [*BASE_BINS, "test.hic.lk.bin", "test.hic.tlb.bin"])
    assert result.hap2_contigs is not None and result.hap2_contigs.name == "test.hic.hap2.p_ctg.gfa"
    await check(
        result,
        {
            "raw_unitigs": ("test.hic.r_utg.gfa", "6a200e4bcb5d8d8fa6ed75874d26edc6"),
            "processed_unitigs": ("test.hic.p_utg.gfa", "6a200e4bcb5d8d8fa6ed75874d26edc6"),
            "primary_contigs": ("test.hic.p_ctg.gfa", "6d0ec1ec65b1d2b295ac94b6fef25a27"),
            "alternate_contigs": None,
            "hap1_contigs": ("test.hic.hap1.p_ctg.gfa", "f67a8fdfa756961360732c79d189054d"),
        },
    )


async def expect_rejected(label: str, **inputs: File) -> None:
    try:
        await hifiasm([await fixture(HIFI)], prefix="test", args="-f0", **inputs)
    except ValueError:
        return
    raise AssertionError(f"hifiasm accepted {label}")


@env.task
async def test_hifiasm_rejects_bad_phasing() -> None:
    yak = await fixture(YAK + "test.yak")
    yak2 = await fixture(YAK + "test2.yak")
    hic1 = await fixture(HIC + "test_1.fastq.gz")
    hic2 = await fixture(HIC + "test_2.fastq.gz")
    await expect_rejected(
        "trio binning with Hi-C",
        paternal_kmer_dump=yak,
        maternal_kmer_dump=yak2,
        hic_read1=hic1,
        hic_read2=hic2,
    )
    await expect_rejected("one k-mer dump", paternal_kmer_dump=yak)
    await expect_rejected("one Hi-C mate", hic_read2=hic2)


tests = [
    test_hifiasm_hifi,
    test_hifiasm_hifi_bins,
    test_hifiasm_hifi_two_files,
    test_hifiasm_hifi_trio,
    test_hifiasm_hifi_hic,
    test_hifiasm_rejects_bad_phasing,
]
