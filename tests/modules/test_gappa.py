"""Tests for flyte_bio.modules.gappa.

The cases are upstream's gappa/examineassign, gappa/examinegraft,
gappa/examineheattree and gappa/editmerge tests (modules commit 9a1cc9f,
gappa 0.9.0), all on upstream's EPA-ng placements of three queries.

- examine assign: with upstream's config's ``--per-query-results --krona
  --sativa``, every output by md5; once more with the taxon file gzipped
  (built in the test, as upstream builds it), to the same md5s.
- examine graft: the grafted tree, by md5.
- examine heat-tree: with upstream's config writing all four tree formats.
  Upstream marks the SVG, colour scale and log unstable, and snapshots the
  SVG's lines 4-8 (its colour gradient) and the colour scale's lines 4-6;
  the other trees by md5.
- edit merge: upstream derives a second jplace from the first by renaming its
  queries (``renamed_<name>``), here in Python the same way, and snapshots
  the merged file's tree, fields and sorted query names, not its md5 (a
  jplace records when it was written). The tree is checked by the md5 of
  upstream's snapshotted tree string.
"""

import gzip
import hashlib
import json
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.gappa import gappa_edit_merge, gappa_examine_assign, gappa_examine_graft, gappa_examine_heat_tree
from tests.framework import assert_md5, env, fixture

JPLACE = "delete_me/gappa/epa_result.jplace.gz"
TAXONOMY = "delete_me/gappa/gappa_taxonomy.tsv"
ASSIGN_ARGS = "--per-query-results --krona --sativa"
ASSIGN_MD5S = {
    "profile": "ce604fa03604462358b69cfbff7c1593",
    "labelled_tree": "39de8a46303ea757f6eb94478db1ffe0",
    "per_query": "16e1fbcc5f79588a0ed21e495cdd91ee",
    "krona": "19110ac998eaa7532f2a92772cd35aa6",
    "sativa": "28020fc2f8440fabb5ae83b0c7114cdb",
}


async def file_bytes(file: File) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read())


def local_file(name: str, data: bytes) -> str:
    path = Path(tempfile.mkdtemp()) / name
    path.write_bytes(data)
    return str(path)


async def check_assign(taxonomy: File) -> None:
    result = await gappa_examine_assign(await fixture(JPLACE), taxonomy, prefix="test", args=ASSIGN_ARGS)
    names = {
        "profile": "test.profile.tsv",
        "labelled_tree": "test.labelled_tree.newick",
        "per_query": "test.per_query.tsv",
        "krona": "test.krona.profile",
        "sativa": "test.sativa.tsv",
    }
    for key, md5 in ASSIGN_MD5S.items():
        out = getattr(result, key)
        assert out is not None, f"gappa wrote no {names[key]}"
        await assert_md5(out, md5, label=names[key])


@env.task
async def test_gappa_examine_assign() -> None:
    await check_assign(await fixture(TAXONOMY))


@env.task
async def test_gappa_examine_assign_gzipped_taxon_file() -> None:
    taxonomy = await file_bytes(await fixture(TAXONOMY))
    await check_assign(await File.from_local(local_file("gappa_taxonomy.tsv.gz", gzip.compress(taxonomy))))


@env.task
async def test_gappa_examine_graft() -> None:
    tree = await gappa_examine_graft(await fixture(JPLACE), prefix="test")
    assert tree.name == "test.epa_result.newick", tree.name
    await assert_md5(tree, "97e54b6a9576a78b774fd63a050e5474", label="test.epa_result.newick")


SVG_LINES = [
    "    <defs>",
    '        <linearGradient id="PaletteGradient_1804289383" x1="0%" y1="100%" x2="0%" y2="0%" spreadMethod="pad" >',
    '            <stop offset="0%" stop-color="#81bfff" stop-opacity="1" />',
    '            <stop offset="50%" stop-color="#c040be" stop-opacity="1" />',
    '            <stop offset="100%" stop-color="#000000" stop-opacity="1" />',
]
COLOUR_LINES = [
    "    At 0.600: Label '0.6', Color #9a3398",
    "    At 0.800: Label '0.8', Color #4d1a4c",
    "    At 1.000: Label '1', Color #000000",
]


@env.task
async def test_gappa_examine_heat_tree() -> None:
    result = await gappa_examine_heat_tree(
        await fixture(JPLACE),
        prefix="test",
        args="--write-newick-tree --write-nexus-tree --write-phyloxml-tree --write-svg-tree",
    )
    for out, name, md5 in [
        (result.newick, "test.tree.newick", "f055e11c9b44b6d3afb502aba99d9579"),
        (result.nexus, "test.tree.nexus", "e1cfdcb63f3a751677d0e78fb6b9230a"),
        (result.phyloxml, "test.tree.phyloxml", "ea35b07e4fd8d174fd6b2f204629cc4e"),
    ]:
        assert out is not None and out.name == name, f"gappa wrote no {name}"
        await assert_md5(out, md5, label=name)
    assert result.svg is not None and result.svg.name == "test.tree.svg", "gappa wrote no test.tree.svg"
    assert result.colours.name == "test.colours.txt" and result.log.name == "test.log"
    svg = (await file_bytes(result.svg)).decode().splitlines()
    assert svg[3:8] == SVG_LINES, svg[3:8]
    colours = (await file_bytes(result.colours)).decode().splitlines()
    assert colours[3:6] == COLOUR_LINES, colours[3:6]


MERGED_FIELDS = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]
MERGED_QUERIES = [
    "A0A1B7XJK0",
    "A0A1J4NSI0",
    "B4F2Z1",
    "renamed_A0A1B7XJK0",
    "renamed_A0A1J4NSI0",
    "renamed_B4F2Z1",
]
# The md5 of upstream's snapshotted merged tree (7514 characters).
MERGED_TREE_MD5 = "3b66a6797089ff7ac6b66dacaee02a9f"


@env.task
async def test_gappa_edit_merge() -> None:
    source = await fixture(JPLACE)
    text = gzip.decompress(await file_bytes(source)).decode()
    renamed = text.replace('"n": ["', '"n": ["renamed_')
    second = await File.from_local(local_file("renamed.jplace.gz", gzip.compress(renamed.encode())))
    out = await gappa_edit_merge([source, second], prefix="test")
    assert out.name == "test.merge.jplace", out.name
    merged = json.loads(await file_bytes(out))
    assert len(merged["placements"]) == 6, len(merged["placements"])
    assert hashlib.md5(merged["tree"].encode()).hexdigest() == MERGED_TREE_MD5, merged["tree"][:80]
    assert merged["fields"] == MERGED_FIELDS, merged["fields"]
    queries = sorted(name for placement in merged["placements"] for name in placement["n"])
    assert queries == MERGED_QUERIES, queries


tests = [
    test_gappa_examine_assign,
    test_gappa_examine_assign_gzipped_taxon_file,
    test_gappa_examine_graft,
    test_gappa_examine_heat_tree,
    test_gappa_edit_merge,
]
