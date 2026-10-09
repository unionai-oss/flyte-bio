"""Tests for flyte_bio.modules.stringtie.

Expected md5s are the upstream stringtie/stringtie snapshot values (no extra
args, as upstream's test config). The transcripts GTF header records the
command line, including ``-p <threads>``.
"""

from flyte_bio.modules.stringtie import StringTieResult, stringtie
from tests.framework import assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"
GTF = "genomics/sarscov2/genome/genome.gtf"
THREADS = 2
BALLGOWN_SAME = {
    "e2t.ctab": "e981c0038295ae54b63cedb1083f1540",
    "i2t.ctab": "8a117c8aa4334b4c2d4711932b006fb4",
    "i_data.ctab": "be3abe09740603213f83d50dcf81427f",
    "t_data.ctab": "3b66c065da73ae0dd41cc332eff6a818",
}


async def run(strandedness: str, with_gtf: bool) -> StringTieResult:
    bam = await fixture(BAM)
    gtf = await fixture(GTF) if with_gtf else None
    return await stringtie(bam, prefix="test", gtf=gtf, strandedness=strandedness, threads=THREADS)


async def check(result: StringTieResult, transcripts: str, abundance: str, ballgown: dict[str, str] | None) -> None:
    await assert_md5(result.transcripts_gtf, transcripts, label="transcripts.gtf")
    await assert_md5(result.abundance, abundance, label="gene.abundance.txt")
    if ballgown:
        tables = {f.path.rsplit("/", 1)[-1]: f async for f in result.results.walk() if "/test.ballgown/" in f.path}
        for name, md5 in ballgown.items():
            assert name in tables, f"missing ballgown {name} (have {sorted(tables)})"
            await assert_md5(tables[name], md5, label=f"ballgown {name}")


@env.task
async def test_stringtie_forward() -> None:
    # upstream case: stringtie/stringtie "sarscov2 [bam] - forward strandedness"
    result = await run("forward", with_gtf=False)
    await check(result, "6087dfc9700a52d9e4a1ae3fcd1d1dfd", "d6f5c8cadb8458f1df0427cf790246e3", None)


@env.task
async def test_stringtie_reverse() -> None:
    # upstream case: stringtie/stringtie "sarscov2 [bam] - reverse strandedness"
    result = await run("reverse", with_gtf=False)
    await check(result, "01d6da00a3c458420841e57427297183", "d6f5c8cadb8458f1df0427cf790246e3", None)


@env.task
async def test_stringtie_forward_reference() -> None:
    # upstream case: stringtie/stringtie "sarscov2 [bam] - forward strandedness + reference annotation"
    result = await run("forward", with_gtf=True)
    await check(
        result,
        "37154e7bda96544f24506ee902bb561d",
        "7d8bce7f2a922e367cedccae7267c22e",
        {**BALLGOWN_SAME, "e_data.ctab": "6b4cf69bc03f3f69890f972a0e8b7471"},
    )


@env.task
async def test_stringtie_reverse_reference() -> None:
    # upstream case: stringtie/stringtie "sarscov2 [bam] - reverse strandedness + reference annotation"
    result = await run("reverse", with_gtf=True)
    await check(
        result,
        "fbabb4e3888bbede67f11f692e484880",
        "7385b870b955dae2c2ab78a70cf05cce",
        {**BALLGOWN_SAME, "e_data.ctab": "879b6696029d19c4737b562e9d149218"},
    )


tests = [
    test_stringtie_forward,
    test_stringtie_reverse,
    test_stringtie_forward_reference,
    test_stringtie_reverse_reference,
]
