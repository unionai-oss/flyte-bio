"""Tests for flyte_bio.modules.odgi.

The cases are upstream's odgi/build, odgi/sort, odgi/unchop, odgi/squeeze,
odgi/view, odgi/layout, odgi/draw, odgi/viz and odgi/stats tests (modules
commit 9a1cc9f, odgi 0.9.0), with odgi/layout's test config's ``ext.args``
(``--out test.lay --tsv test.tsv``, here ``tsv=True``). Upstream marks every
graph, layout, image and stats output unstable (odgi's sort and layout are
stochastic and multithreaded) and checks only their names; these tests also
check that graphs and layouts are nonempty and images are PNGs. odgi stats
without options prints nothing (upstream notes its output is empty), so that
case checks the name only.

Upstream tests each tool on fixtures; :func:`test_odgi_pangenome_chain`
also builds the graph from the pangenome FASTA, as a pangenome pipeline does:
wfmash → seqwish → smoothxg → gfaffix → odgi build → sort → unchop →
layout → draw, with viz, stats, view and squeeze of the results, checking
that the paths (one per input sequence, plus the consensus paths smoothxg
adds) survive every step.
"""

from flyte.io import File

from flyte_bio.modules.gfaffix import gfaffix
from flyte_bio.modules.odgi import (
    odgi_build,
    odgi_draw,
    odgi_layout,
    odgi_sort,
    odgi_squeeze,
    odgi_stats,
    odgi_unchop,
    odgi_view,
    odgi_viz,
)
from flyte_bio.modules.seqwish import seqwish_induce
from flyte_bio.modules.smoothxg import smoothxg
from flyte_bio.modules.wfmash import wfmash
from tests.framework import assert_nonempty, env, fixture

PANGENOME = "pangenomics/homo_sapiens/"
SMOOTHXG_GFA = PANGENOME + "pangenome.smoothxg.gfa"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


async def odgi_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


def check_name(file: File, name: str) -> None:
    assert file.path.endswith("/" + name), file.path


async def check_png(file: File, name: str) -> None:
    check_name(file, name)
    async with file.open("rb") as fh:
        head = bytes(await fh.read(len(PNG_SIGNATURE)))
    assert head == PNG_SIGNATURE, f"{name}: starts with {head!r}"


async def check_graph(file: File, name: str) -> None:
    check_name(file, name)
    await assert_nonempty(file, label=name)


async def gfa_path_names(gfa: File) -> list[str]:
    return sorted(line.split("\t")[1] for line in (await odgi_text(gfa)).splitlines() if line.startswith("P\t"))


def input_paths(paths: list[str]) -> list[str]:
    """``paths`` without the consensus paths smoothxg adds (``Consensus_<n>``)."""
    return [p for p in paths if not p.startswith("Consensus_")]


@env.task
async def test_odgi_build() -> None:
    # upstream case: odgi/build "homo_sapiens - pangenome - pangenome_smoothxg_gfa"
    await check_graph(await odgi_build(await fixture(SMOOTHXG_GFA), prefix="test"), "test.og")


@env.task
async def test_odgi_sort() -> None:
    # upstream case: odgi/sort "homo_sapiens - pangenome - pangenome_smoothxg_gfa"
    await check_graph(await odgi_sort(await fixture(SMOOTHXG_GFA), prefix="test"), "test.og")


@env.task
async def test_odgi_unchop() -> None:
    # upstream case: odgi/unchop "homo_sapiens - pangenome - pangenome_smoothxg_gfa"
    await check_graph(await odgi_unchop(await fixture(SMOOTHXG_GFA), prefix="test"), "test.og")


@env.task
async def test_odgi_squeeze() -> None:
    # upstream case: odgi/squeeze "homo_sapiens - pangenome - odgi - pangenome_og"
    graph = await fixture(PANGENOME + "odgi/pangenome.og")
    await check_graph(await odgi_squeeze([graph], prefix="test"), "test.og")


@env.task
async def test_odgi_view() -> None:
    # upstream case: odgi/view "homo_sapiens - pangenome - pangenome_smoothxg_gfa"
    gfa = await fixture(SMOOTHXG_GFA)
    view = await odgi_view(gfa, prefix="test")
    check_name(view, "test.gfa")
    assert await gfa_path_names(view) == await gfa_path_names(gfa), "odgi view changed the paths"


@env.task
async def test_odgi_layout() -> None:
    # upstream case: odgi/layout "homo_sapiens - pangenome - pangenome_smoothxg_gfa"
    result = await odgi_layout(await fixture(SMOOTHXG_GFA), prefix="test", tsv=True)
    assert result.tsv is not None
    await check_graph(result.lay, "test.lay")
    check_name(result.tsv, "test.tsv")
    header = (await odgi_text(result.tsv)).splitlines()[0].split("\t")
    assert header[:3] == ["idx", "X", "Y"], header


@env.task
async def test_odgi_draw() -> None:
    # upstream case: odgi/draw "homo_sapiens - pangenome - odgi - pangenome_og - pangenome_lay"
    graph = await fixture(PANGENOME + "odgi/pangenome.og")
    layout = await fixture(PANGENOME + "odgi/pangenome.lay")
    await check_png(await odgi_draw(graph, layout, prefix="test"), "test.png")


@env.task
async def test_odgi_viz() -> None:
    # upstream case: odgi/viz "homo_sapiens - pangenome - pangenome_smoothxg_gfa"
    await check_png(await odgi_viz(await fixture(SMOOTHXG_GFA), prefix="test"), "test.png")


@env.task
async def test_odgi_stats() -> None:
    # upstream case: odgi/stats "homo_sapiens - pangenome - pangenome_smoothxg_gfa"
    check_name(await odgi_stats(await fixture(SMOOTHXG_GFA), prefix="test"), "test.og.stats.tsv")


@env.task
async def test_odgi_pangenome_chain() -> None:
    fasta_gz = await fixture(PANGENOME + "pangenome.fa.gz")
    gzi = await fixture(PANGENOME + "pangenome.fa.gz.gzi")
    fai = await fixture(PANGENOME + "pangenome.fa.gz.fai")
    fasta = await fixture(PANGENOME + "pangenome.fa")
    sequences = sorted(line.split("\t")[0] for line in (await odgi_text(fai)).splitlines())

    mappings = await wfmash(fasta_gz, gzi, fai, query_self=True, prefix="chain", args="-s 500 -X")
    induced = await seqwish_induce([mappings.paf], fasta, prefix="chain")
    smoothed = await smoothxg(induced.gfa, prefix="chain", args="--n-haps 10")
    normalised = await gfaffix(smoothed.gfa, prefix="chain")
    assert await gfa_path_names(induced.gfa) == sequences, "seqwish: paths differ from the input sequences"
    for label, gfa in (("smoothxg", smoothed.gfa), ("gfaffix", normalised.gfa)):
        paths = await gfa_path_names(gfa)
        assert input_paths(paths) == sequences, f"{label}: paths {paths} differ from the input sequences"
        assert len(paths) > len(sequences), f"{label}: no consensus paths"

    built = await odgi_build(normalised.gfa, prefix="chain")
    ordered = await odgi_sort(built, prefix="chain.sorted", args="-p Ygs")
    unchopped = await odgi_unchop(ordered, prefix="chain.unchopped")
    layout = await odgi_layout(unchopped, prefix="chain", tsv=True)
    await check_graph(layout.lay, "chain.lay")
    await check_png(await odgi_draw(unchopped, layout.lay, prefix="chain.draw"), "chain.draw.png")
    await check_png(await odgi_viz(unchopped, prefix="chain.viz"), "chain.viz.png")

    view = await odgi_view(unchopped, prefix="chain")
    view_paths = await gfa_path_names(view)
    assert view_paths == await gfa_path_names(normalised.gfa), "odgi build to view: paths changed"
    squeezed = await odgi_squeeze([built, ordered], prefix="chain.squeezed")
    squeezed_view = await odgi_view(squeezed, prefix="chain.squeezed")
    squeezed_paths = await gfa_path_names(squeezed_view)
    assert set(squeezed_paths) == set(view_paths), "odgi squeeze: paths differ from the inputs"

    stats = await odgi_stats(unchopped, prefix="chain", args="-S")
    check_name(stats, "chain.og.stats.tsv")
    rows = [line.split("\t") for line in (await odgi_text(stats)).splitlines()]
    assert rows[0][:5] == ["#length", "nodes", "edges", "paths", "steps"], rows[0]
    assert int(rows[1][3]) == len(view_paths), rows


tests = [
    test_odgi_build,
    test_odgi_sort,
    test_odgi_unchop,
    test_odgi_squeeze,
    test_odgi_view,
    test_odgi_layout,
    test_odgi_draw,
    test_odgi_viz,
    test_odgi_stats,
    test_odgi_pangenome_chain,
]
