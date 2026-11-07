"""Tests for flyte_bio.modules.riker.

The cases and expected md5s are upstream's riker/multi tests (modules commit
9a1cc9f, riker 0.5.0), each with its ``--tools`` from upstream's test
parameters. As upstream, charts (and in the RNA case the biotype and insert
size histogram tables) are checked by name only, and every case checks the
full set of files written. The module's 4 CPUs match upstream's test CPUs.

Upstream's two hybcap failure cases fail because baits come without targets
(or targets without baits); here :func:`riker_multi` refuses those before
running.
"""

from flyte.io import Dir, File

from flyte_bio.modules.riker import riker_multi
from tests.framework import assert_md5, env, fixture

SARSCOV2_BAM = "genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"
SARSCOV2_FASTA = "genomics/sarscov2/genome/genome.fasta"
BAITS = "genomics/sarscov2/genome/picard/baits.interval_list"
TARGETS = "genomics/sarscov2/genome/picard/targets.interval_list"

ALIGNMENT = {
    "test.alignment-metrics.txt": "ef2421780f4e1127dac0ae3d6702d33f",
    "test.base-distribution-by-cycle.txt": "d6e323db9bbabde773842be47db63ccc",
    "test.mean-quality-by-cycle.txt": "476932952a367fe567eba75d6ff51630",
    "test.quality-score-distribution.txt": "d5f17682727e31f05dc8b29b7c06b3ab",
    "test.base-distribution-by-cycle.pdf": None,
    "test.mean-quality-by-cycle.pdf": None,
    "test.quality-score-distribution.pdf": None,
}
ISIZE = {
    "test.isize-histogram.txt": "68742e205c83b880802c3c5bbe7e882c",
    "test.isize-metrics.txt": "c3bef6653b685c99e373a1dfccedc52d",
    "test.isize-histogram.pdf": None,
}
WGS_GCBIAS = {
    "test.gcbias-detail.txt": "44cd58d1a0bd81ca8ea5d3fa7c07e4db",
    "test.gcbias-summary.txt": "4915d05211aaae398fc0f29ab3cc79b8",
    "test.wgs-coverage.txt": "8fd648b0045ff06fc3ba46772e33255c",
    "test.wgs-metrics.txt": "d7bb5636d4abeefe85f5ca2a6674adb2",
    "test.gcbias-chart.pdf": None,
    "test.wgs-coverage.pdf": None,
}
ERROR = {
    "test.error-indel.txt": "9626acfbfc0c901c8a2999df3796cae0",
    "test.error-mismatch.txt": "e5a53c4bd345321a872e94752eabae18",
    "test.error-overlap.txt": "4426d9733ae3e59dc2496b0c247c22f2",
}
HYBCAP = {"test.hybcap-metrics.txt": "65d0a112bd3508d6c3b5328c8fcfcbdb"}
ERROR_ARGS = "--error::stratify-by bq isize gc"


async def check(results: Dir, expected: dict[str, str | None]) -> None:
    """``results`` holds exactly the ``expected`` files; those with an md5 match it."""
    files = {f.name: f for f in await results.list_files()}
    assert set(files) == set(expected), f"missing {set(expected) - set(files)}, extra {set(files) - set(expected)}"
    for name, md5 in expected.items():
        if md5 is not None:
            await assert_md5(files[name], md5, label=name)


async def sarscov2_bam() -> tuple[File, File]:
    return await fixture(SARSCOV2_BAM), await fixture(SARSCOV2_BAM + ".bai")


async def sarscov2_fasta() -> tuple[File, File]:
    return await fixture(SARSCOV2_FASTA), await fixture(SARSCOV2_FASTA + ".fai")


@env.task
async def test_riker_alignment_basic_isize() -> None:
    # upstream case: "sarscov2 - paired_end - bam - alignment basic isize"
    bam, bai = await sarscov2_bam()
    results = await riker_multi(bam, "test", alignment_index=bai, args="--tools alignment basic isize")
    await check(results, ALIGNMENT | ISIZE)


@env.task
async def test_riker_wgs_gcbias() -> None:
    # upstream case: "sarscov2 - paired_end - bam - wgs gcbias"
    (bam, bai), (fasta, fai) = await sarscov2_bam(), await sarscov2_fasta()
    results = await riker_multi(
        bam, "test", alignment_index=bai, fasta=fasta, fai=fai, args="--tools wgs gcbias --wgs::per-base-coverage"
    )
    await check(results, WGS_GCBIAS | {"test.wgs-per-base.bw": "3ccd00aa7a7e120311f50e1a71c1cb8f"})


@env.task
async def test_riker_wgs_gcbias_alignment_basic_isize() -> None:
    # upstream case: "sarscov2 - paired_end - bam - wgs gcbias alignment basic isize"
    (bam, bai), (fasta, fai) = await sarscov2_bam(), await sarscov2_fasta()
    results = await riker_multi(
        bam, "test", alignment_index=bai, fasta=fasta, fai=fai, args="--tools wgs gcbias alignment basic isize"
    )
    await check(results, ALIGNMENT | ISIZE | WGS_GCBIAS)


@env.task
async def test_riker_cram_alignment_basic_isize() -> None:
    # upstream case: "homo_sapiens - paired_end - cram - alignment basic isize"
    cram = await fixture("genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram")
    crai = await fixture("genomics/homo_sapiens/illumina/cram/test.paired_end.sorted.cram.crai")
    fasta = await fixture("genomics/homo_sapiens/genome/genome.fasta")
    fai = await fixture("genomics/homo_sapiens/genome/genome.fasta.fai")
    results = await riker_multi(
        cram, "test", alignment_index=crai, fasta=fasta, fai=fai, args="--tools alignment basic isize"
    )
    expected = dict.fromkeys(ALIGNMENT | ISIZE) | {
        "test.alignment-metrics.txt": "2e8b897e737afd500c17de8526ea055b",
        "test.base-distribution-by-cycle.txt": "5fa28912c8eb847f59f01c22cefaf872",
        "test.isize-histogram.txt": "fc8e1ced54e36439132de087895d85fc",
        "test.isize-metrics.txt": "5abd96d5d46783ff122b110b95b7eedd",
        "test.mean-quality-by-cycle.txt": "70c31ea94a4a6be25a12a1a2221b7ffa",
        "test.quality-score-distribution.txt": "1e4ceaf9de78390f1f542eec5dc7275f",
    }
    await check(results, expected)


@env.task
async def test_riker_hybcap() -> None:
    # upstream case: "sarscov2 - paired_end - bam - hybcap"
    (bam, bai), (fasta, fai) = await sarscov2_bam(), await sarscov2_fasta()
    baits, targets = await fixture(BAITS), await fixture(TARGETS)
    results = await riker_multi(
        bam,
        "test",
        alignment_index=bai,
        fasta=fasta,
        fai=fai,
        hybcap_baits=baits,
        hybcap_targets=targets,
        args="--tools hybcap",
    )
    await check(results, HYBCAP)


@env.task
async def test_riker_error() -> None:
    # upstream case: "sarscov2 - paired_end - bam - error"
    (bam, bai), (fasta, fai) = await sarscov2_bam(), await sarscov2_fasta()
    results = await riker_multi(
        bam, "test", alignment_index=bai, fasta=fasta, fai=fai, args=f"--tools error {ERROR_ARGS}"
    )
    await check(results, ERROR)


@env.task
async def test_riker_error_alignment() -> None:
    # upstream case: "sarscov2 - paired_end - bam - error alignment"
    (bam, bai), (fasta, fai) = await sarscov2_bam(), await sarscov2_fasta()
    results = await riker_multi(
        bam, "test", alignment_index=bai, fasta=fasta, fai=fai, args=f"--tools error alignment {ERROR_ARGS}"
    )
    await check(results, ERROR | {"test.alignment-metrics.txt": ALIGNMENT["test.alignment-metrics.txt"]})


@env.task
async def test_riker_all_tools() -> None:
    # upstream case: "sarscov2 - paired_end - bam - all tools"
    (bam, bai), (fasta, fai) = await sarscov2_bam(), await sarscov2_fasta()
    baits, targets = await fixture(BAITS), await fixture(TARGETS)
    results = await riker_multi(
        bam,
        "test",
        alignment_index=bai,
        fasta=fasta,
        fai=fai,
        hybcap_baits=baits,
        hybcap_targets=targets,
        args=f"--tools alignment basic isize wgs gcbias hybcap error {ERROR_ARGS}",
    )
    await check(results, ALIGNMENT | ISIZE | WGS_GCBIAS | ERROR | HYBCAP)


async def refused(**hybcap: File) -> bool:
    bam, bai = await sarscov2_bam()
    try:
        await riker_multi(bam, "test", alignment_index=bai, args="--tools hybcap", **hybcap)
    except ValueError:
        return True
    return False


@env.task
async def test_riker_baits_without_targets_fails() -> None:
    # upstream case: "sarscov2 - paired_end - bam - baits without targets fails"
    assert await refused(hybcap_baits=await fixture(BAITS)), "riker_multi ran with baits but no targets"


@env.task
async def test_riker_targets_without_baits_fails() -> None:
    # upstream case: "sarscov2 - paired_end - bam - targets without baits fails"
    assert await refused(hybcap_targets=await fixture(TARGETS)), "riker_multi ran with targets but no baits"


@env.task
async def test_riker_rna() -> None:
    # upstream case: "homo_sapiens - paired_end - bam - rna"
    bam = await fixture("genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam")
    bai = await fixture("genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam.bai")
    gff = await fixture("genomics/homo_sapiens/genome/genome.gff3")
    results = await riker_multi(bam, "test", alignment_index=bai, rna_gene_model=gff, args="--tools rna")
    await check(
        results,
        {
            "test.rna-insert-size.txt": "1ba71842d0e1f52b9a92067126ed5ec1",
            "test.rna-metrics.txt": "bf9772f8a31ded0f3e7e3102204394cd",
            "test.rna-biotype.txt": None,
            "test.rna-insert-size-histogram.txt": None,
            "test.rna-coverage.pdf": None,
        },
    )


tests = [
    test_riker_alignment_basic_isize,
    test_riker_wgs_gcbias,
    test_riker_wgs_gcbias_alignment_basic_isize,
    test_riker_cram_alignment_basic_isize,
    test_riker_hybcap,
    test_riker_error,
    test_riker_error_alignment,
    test_riker_all_tools,
    test_riker_baits_without_targets_fails,
    test_riker_targets_without_baits_fails,
    test_riker_rna,
]
