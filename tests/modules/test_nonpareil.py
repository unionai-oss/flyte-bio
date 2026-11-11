"""Tests for flyte_bio.modules.nonpareil.

The cases are upstream's nonpareil/nonpareil, nonpareil/curve, nonpareil/set
and nonpareil/nonpareilcurvesr tests (modules commit 9a1cc9f, Nonpareil
3.5.5), with their configs' options (``-X 100``; the plots' colours). The
reads are Candidatus Portiera aleyrodidarum's, gzipped or decompressed with
the gunzip module as upstream's setup. As upstream (the outputs carry
run-dependent values), the outputs are checked by name, the ``.npo``'s first
line, the log's verdict, the summaries' contents, and the plots' presence.
"""

from flyte.io import File

from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.nonpareil import (
    NonpareilResult,
    nonpareil_curve,
    nonpareil_nonpareil,
    nonpareil_nonpareilcurvesr,
    nonpareil_set,
)
from tests.framework import env, fixture

READS = "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/illumina/fastq/test_1.fastq.gz"


async def nonpareil_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check_nonpareil(result: NonpareilResult) -> None:
    assert (result.npa.name, result.npc.name) == ("test.npa", "test.npc"), (result.npa.name, result.npc.name)
    first = (await nonpareil_text(result.npo)).splitlines()[0]
    assert first == "# @impl: Nonpareil", f"test.npo's first line: {first!r}"
    last = (await nonpareil_text(result.npl)).splitlines()[-1]
    assert "Everything seems correct" in last, f"test.npl's last line: {last!r}"


@env.task
async def test_nonpareil_nonpareil_gzipped() -> None:
    # upstream case: nonpareil/nonpareil "candidatus_portiera_aleyrodidarum - gzipped"
    await check_nonpareil(
        await nonpareil_nonpareil(await fixture(READS), "fastq", "kmer", prefix="test", args="-X 100")
    )


@env.task
async def test_nonpareil_nonpareil_uncompressed() -> None:
    # upstream case: nonpareil/nonpareil "candidatus_portiera_aleyrodidarum - uncompressed"
    reads = await gunzip(await fixture(READS))
    await check_nonpareil(await nonpareil_nonpareil(reads, "fastq", "kmer", prefix="test", args="-X 100"))


@env.task
async def test_nonpareil_curve() -> None:
    # upstream case: nonpareil/curve "candidatus_portiera_aleyrodidarum"
    result = await nonpareil_nonpareil(await fixture(READS), "fastq", "kmer", prefix="test", args="-X 100")
    png = await nonpareil_curve(result.npo, prefix="test", r_args="col = 'green'")
    assert png.name == "test.png", png.name


@env.task
async def test_nonpareil_set() -> None:
    # upstream case: nonpareil/set "candidatus_portiera_aleyrodidarum"
    reads = await fixture(READS)
    npos = [(await nonpareil_nonpareil(reads, "fastq", "kmer", prefix=p, args="-X 100")).npo for p in ("test", "test2")]
    png = await nonpareil_set(npos, prefix="test", r_args="col = c('blue', 'red')")
    assert png.name == "test.png", png.name


@env.task
async def test_nonpareil_nonpareilcurvesr() -> None:
    # upstream case: nonpareil/nonpareilcurvesr "candidatus_portiera_aleyrodidarum"
    reads = await fixture(READS)
    npos = [(await nonpareil_nonpareil(reads, "fastq", "kmer", prefix=p)).npo for p in ("test", "test2")]
    result = await nonpareil_nonpareilcurvesr(npos, prefix="test")
    assert 'modelR": 0.' in await nonpareil_text(result.json), "test.json has no modelR"
    assert "test\t0." in await nonpareil_text(result.tsv), "test.tsv has no row for test"
    assert 'test",0.' in await nonpareil_text(result.csv), "test.csv has no row for test"
    assert result.pdf.name == "test.pdf", result.pdf.name


tests = [
    test_nonpareil_nonpareil_gzipped,
    test_nonpareil_nonpareil_uncompressed,
    test_nonpareil_curve,
    test_nonpareil_set,
    test_nonpareil_nonpareilcurvesr,
]
