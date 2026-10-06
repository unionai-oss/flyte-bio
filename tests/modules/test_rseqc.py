"""Tests for flyte_bio.modules.rseqc.

One test per upstream rseqc/* module case (sarscov2 paired-end sorted BAM).
Expected md5s are the upstream snapshot values; PDFs aren't byte-stable, so
— as upstream — only their names are checked.
"""

from flyte.io import Dir

from flyte_bio.modules import rseqc
from tests.framework import assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"
BED = "genomics/sarscov2/genome/bed/test.bed"
BED12 = "genomics/sarscov2/genome/bed/test.bed12"


async def inputs(bed: str | None = None):
    bam, bai = await fixture(BAM), await fixture(BAM + ".bai")
    return bam, bai, (await fixture(bed) if bed else None)


async def check(results: Dir, md5s: dict[str, str], present: tuple[str, ...] = ()) -> None:
    for name, md5 in md5s.items():
        f = await results.get_file(name)
        assert f is not None, f"missing {name}"
        await assert_md5(f, md5, label=name)
    for name in present:
        assert await results.get_file(name) is not None, f"missing {name}"


@env.task
async def test_rseqc_bam_stat() -> None:
    # upstream case: rseqc/bamstat "sarscov2 - [meta] - bam"
    bam, bai, _ = await inputs()
    out = await rseqc.bam_stat(bam=bam, bai=bai, prefix="test")
    await check(out, {"test.bam_stat.txt": "2675857864c1d1139b2a19d25dc36b09"})


@env.task
async def test_rseqc_infer_experiment() -> None:
    # upstream case: rseqc/inferexperiment "sarscov2 - [[meta] - bam] - bed"
    bam, bai, bed = await inputs(BED)
    out = await rseqc.infer_experiment(bam=bam, bai=bai, bed=bed, prefix="test")
    await check(out, {"test.infer_experiment.txt": "f9d0bfc239df637cd8aeda40ade3c59a"})


@env.task
async def test_rseqc_inner_distance() -> None:
    # upstream case: rseqc/innerdistance "sarscov2 - [[meta] - bam] - bed"
    bam, bai, bed = await inputs(BED12)
    out = await rseqc.inner_distance(bam=bam, bai=bai, bed=bed, prefix="test", single_end=False)
    await check(
        out,
        {
            "test.inner_distance.txt": "a1acc9def0f64a5500d4c4cb47cbe32b",
            "test.inner_distance_freq.txt": "3fc037501f5899b5da009c8ce02fc25e",
            "test.inner_distance_mean.txt": "58398b7d5a29a5e564f9e3c50b55996c",
            "test.inner_distance_plot.r": "5859fbd5b42046d47e8b9aa85077f4ea",
        },
        present=("test.inner_distance_plot.pdf",),
    )


@env.task
async def test_rseqc_junction_annotation() -> None:
    # upstream case: rseqc/junctionannotation "sarscov2 - paired end [bam]"
    bam, bai, bed = await inputs(BED12)
    out = await rseqc.junction_annotation(bam=bam, bai=bai, bed=bed, prefix="test")
    await check(out, {"test.junction_annotation.log": "d75e0f5d62fada8aa9449991b209554c"})


@env.task
async def test_rseqc_junction_saturation() -> None:
    # upstream case: rseqc/junctionsaturation "sarscov2 paired-end [bam]"
    bam, bai, bed = await inputs(BED12)
    out = await rseqc.junction_saturation(bam=bam, bai=bai, bed=bed, prefix="test")
    await check(
        out,
        {"test.junctionSaturation_plot.r": "caa6e63dcb477aabb169882b2f30dadd"},
        present=("test.junctionSaturation_plot.pdf",),
    )


@env.task
async def test_rseqc_read_distribution() -> None:
    # upstream case: rseqc/readdistribution "sarscov2 paired-end [bam]"
    bam, bai, bed = await inputs(BED12)
    out = await rseqc.read_distribution(bam=bam, bai=bai, bed=bed, prefix="test")
    await check(out, {"test.read_distribution.txt": "56893fdc0809d968629a363551a1655f"})


@env.task
async def test_rseqc_read_duplication() -> None:
    # upstream case: rseqc/readduplication "sarscov2 paired-end [bam]"
    bam, bai, _ = await inputs()
    out = await rseqc.read_duplication(bam=bam, bai=bai, prefix="test")
    await check(
        out,
        {
            "test.pos.DupRate.xls": "a859bc2031d46bf1cc4336205847caa3",
            "test.DupRate_plot.r": "3c0325095cee4835b921e57d61c23dca",
            "test.seq.DupRate.xls": "ee8783399eec5a18522a6f08bece338b",
        },
        present=("test.DupRate_plot.pdf",),
    )


tests = [
    test_rseqc_bam_stat,
    test_rseqc_infer_experiment,
    test_rseqc_inner_distance,
    test_rseqc_junction_annotation,
    test_rseqc_junction_saturation,
    test_rseqc_read_distribution,
    test_rseqc_read_duplication,
]
