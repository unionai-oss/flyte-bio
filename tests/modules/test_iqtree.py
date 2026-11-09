"""Tests for flyte_bio.modules.iqtree.

The cases are upstream's iqtree tests (modules commit 9a1cc9f, IQ-TREE 3.1.3),
with their test configs' args. Most of IQ-TREE's outputs depend on its random
seed, so upstream snapshots only their first lines (or that they name a
taxon), the outputs that aren't written, and the md5s of the site-wise
outputs of the "optional" case, which evaluates a fixed tree; these tests
assert the same. As upstream, the "sup" case passes the tree's copy as
``hydrogenase.newick``.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.iqtree import IqtreeResult, iqtree
from tests.framework import assert_md5, env, fixture

SETOXIN = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin.ref"
PHYLOPLACE = "https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/"
HYDROGENASE_TREE = PHYLOPLACE + "PF14720_seed.ft.LGCAT.newick"
HYDROGENASE_ALIGNMENT = PHYLOPLACE + "PF14720_seed.alnfaa"
REPORT_FIRST_LINE = "IQ-TREE 3.1.3 built Jul 26 2026"
LOG_FIRST_LINE = "IQ-TREE version 3.1.3 for Linux x86 64-bit built Jul 26 2026"
OUTPUTS = [
    "phylogeny",
    "report",
    "mldist",
    "lmap_svg",
    "lmap_eps",
    "lmap_quartetlh",
    "sitefreq_out",
    "bootstrap",
    "state",
    "contree",
    "nex",
    "splits",
    "suptree",
    "alninfo",
    "partlh",
    "siteprob",
    "sitelh",
    "treels",
    "rate",
    "mlrate",
    "exch_matrix",
]


async def lines(file: File | None, label: str) -> list[str]:
    assert file is not None, f"iqtree wrote no {label}"
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


def check_absent(result: IqtreeResult, written: set[str]) -> None:
    """Assert that every output upstream's case expects empty wasn't written."""
    for name in OUTPUTS:
        if name in written:
            continue
        value = getattr(result, name)
        assert value in (None, []), f"iqtree wrote unexpected {name}: {value}"


async def check_first_line(file: File | None, label: str, expected: str) -> None:
    first = (await lines(file, label))[0]
    assert first == expected, f"{label} starts {first!r}, expected {expected!r}"


async def check_contains_first_line(file: File | None, label: str, wanted: str) -> None:
    first = (await lines(file, label))[0]
    assert wanted in first, f"{label}'s first line {first!r} lacks {wanted!r}"


async def check_setoxin(result: IqtreeResult) -> None:
    await check_contains_first_line(result.phylogeny, "treefile", "1apf")
    await check_first_line(result.mldist, "mldist", "5")
    await check_first_line(result.report, "report", REPORT_FIRST_LINE)
    await check_first_line(result.log, "log", LOG_FIRST_LINE)


@env.task
async def test_iqtree_setoxin_basic() -> None:
    result = await iqtree(await fixture(SETOXIN), prefix="test")
    check_absent(result, {"phylogeny", "report", "mldist"})
    await check_setoxin(result)


@env.task
async def test_iqtree_setoxin_treels() -> None:
    result = await iqtree(await fixture(SETOXIN), prefix="test", args="-m LG -wt")
    check_absent(result, {"phylogeny", "report", "mldist", "treels"})
    await check_setoxin(result)
    await check_contains_first_line(result.treels, "treels", "1apf")


@env.task
async def test_iqtree_setoxin_bootstrap() -> None:
    result = await iqtree(await fixture(SETOXIN), prefix="test", args="-bb 1000 -wbt -wsplits -m LG")
    check_absent(result, {"phylogeny", "report", "mldist", "bootstrap", "contree", "splits", "nex"})
    await check_setoxin(result)
    await check_contains_first_line(result.bootstrap, "ufboot", "1apf")
    await check_contains_first_line(result.contree, "contree", "1apf")
    await check_contains_first_line(result.splits, "splits", "*")
    assert len(result.nex) == 1, f"iqtree wrote nex files {[f.name for f in result.nex]}, expected one"
    expected = [
        "#nexus",
        "",
        "BEGIN Taxa;",
        "DIMENSIONS ntax=5;",
        "TAXLABELS",
        "[1] '1apf'",
        "[2] '1ahl'",
        "[3] '1atx'",
        "[4] '1sh1'",
        "[5] '1bds'",
        ";",
        "END; [Taxa]",
        "",
    ]
    head = (await lines(result.nex[0], result.nex[0].name))[:13]
    assert head == expected, f"{result.nex[0].name} starts {head}, expected {expected}"


@env.task
async def test_iqtree_hydrogenase_sup() -> None:
    tree = await fixture(HYDROGENASE_TREE)
    # Upstream collects the tree's text into hydrogenase.newick.
    copy = Path(tempfile.mkdtemp()) / "hydrogenase.newick"
    async with tree.open("rb") as fh:
        copy.write_bytes(bytes(await fh.read()))
    result = await iqtree(tree=tree, suptree=await File.from_local(str(copy)), prefix="test")
    check_absent(result, {"suptree"})
    await check_first_line(result.log, "log", LOG_FIRST_LINE)
    assert result.suptree is not None, "iqtree wrote no suptree"
    await assert_md5(result.suptree, "4ad04b2105448a802d5eff346af67786", label="test.suptree")


@env.task
async def test_iqtree_hydrogenase_optional() -> None:
    result = await iqtree(
        await fixture(HYDROGENASE_ALIGNMENT),
        tree_te=await fixture(HYDROGENASE_TREE),
        prefix="test",
        args="-m C10 -asr -alninfo -wpl -wspr -wsl -wsr --mlrate -lmap 150 -wql --mlrate",
    )
    for name in ["mldist", "sitefreq_out", "bootstrap", "contree", "nex", "suptree", "partlh", "exch_matrix"]:
        assert getattr(result, name) in (None, []), f"iqtree wrote unexpected {name}"
    # Upstream snapshots treels and splits as not written.
    assert result.treels is None and result.splits is None, "iqtree wrote unexpected treels or splits"
    await check_contains_first_line(result.phylogeny, "treefile", "E1QGS5_DESB2/226-309")
    eps = (await lines(result.lmap_eps, "lmap.eps"))[:5]
    assert eps == [
        "%!PS-Adobe-3.0 EPSF-3.0",
        "%%BoundingBox: 60 210 550 650",
        "%%Pages: 1",
        "%%Creator: IQ-TREE/TREE-PUZZLE",
        "%%Title: Likelihood Mapping Analysis",
    ], f"lmap.eps starts {eps}"
    svg = (await lines(result.lmap_svg, "lmap.svg"))[:6]
    assert svg == [
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
        '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">',
        "<svg",
        '   xmlns:svg="http://www.w3.org/2000/svg"',
        '   xmlns="http://www.w3.org/2000/svg"',
        '   xmlns:xlink="http://www.w3.org/1999/xlink"',
    ], f"lmap.svg starts {svg}"
    await check_first_line(
        result.lmap_quartetlh, "lmap.quartetlh", "SeqIDs\tlh1\tlh2\tlh3\tweight1\tweight2\tweight3\tarea\tcorner"
    )
    mlrate = (await lines(result.mlrate, "mlrate"))[:6]
    assert mlrate == [
        "# Site-specific subtitution rates determined by maximum likelihood",
        "# This file can be read in MS Excel or in R with command:",
        "#   tab=read.table('test.mlrate',header=TRUE)",
        "# Columns are tab-separated with following meaning:",
        "#   Site:   Alignment site ID",
        "#   Rate:   Site rate estimated by maximum likelihood",
    ], f"mlrate starts {mlrate}"
    await check_first_line(result.report, "report", REPORT_FIRST_LINE)
    await check_first_line(result.log, "log", LOG_FIRST_LINE)
    for name, md5 in [
        ("alninfo", "3636d762d4cdeefc20f8fd3fe76ba0b7"),
        ("siteprob", "9b3cb9648ca042951cac869bb856c117"),
        ("sitelh", "74b2c736eeede06b2ca4e1f60f3e8dab"),
        ("rate", "a142257cce95dffde8dec4e68995e553"),
        ("state", "fc636a9b1b84b02b2ab631d529b4f1d2"),
    ]:
        file = getattr(result, name)
        assert file is not None, f"iqtree wrote no {name}"
        await assert_md5(file, md5, label=f"test.{name}")


tests = [
    test_iqtree_setoxin_basic,
    test_iqtree_setoxin_treels,
    test_iqtree_setoxin_bootstrap,
    test_iqtree_hydrogenase_sup,
    test_iqtree_hydrogenase_optional,
]
