"""Tests for flyte_bio.modules.csvtk.

The cases and expected md5s are upstream's csvtk/join and csvtk/concat tests
(modules commit 9a1cc9f, csvtk 0.37.0), on the bacass test samplesheets (which
are tab-separated despite their ``.csv`` names). Upstream's join ``ext.args``
``--fields 'ID;ID' -p -e -d "\\t" -D ","`` become ``fields``, ``delimiter`` and
``out_delimiter`` inputs plus ``args="-p -e"``.
"""

from flyte.io import File

from flyte_bio.modules.csvtk import csvtk_concat, csvtk_join
from tests.framework import assert_md5, env, fixture

BACASS = "https://github.com/nf-core/test-datasets/raw/bacass/"


def name_of(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def bacass(name: str) -> File:
    return await fixture(BACASS + name)


@env.task
async def test_csvtk_join_csv() -> None:
    tables = [await bacass("bacass_hybrid.csv"), await bacass("bacass_short.csv")]
    out = await csvtk_join(tables, prefix="test", fields="ID;ID", delimiter="\t", out_delimiter=",", args="-p -e")
    assert name_of(out) == "test.csv", f"expected test.csv, got {name_of(out)}"
    await assert_md5(out, "0ab28b22412f9d4bfb0cf7de262746cb", label="test.csv")


@env.task
async def test_csvtk_join_tsv_output() -> None:
    tables = [await bacass("bacass_hybrid.csv"), await bacass("bacass_short.csv")]
    out = await csvtk_join(tables, prefix="test", fields="ID;ID", delimiter="\t", args="-p -e -T")
    assert name_of(out) == "test.tsv", f"expected test.tsv, got {name_of(out)}"
    await assert_md5(out, "ad4dbe79340e364207698f870ec58cfd", label="test.tsv")


@env.task
async def test_csvtk_concat_tsv_to_csv() -> None:
    short = await bacass("bacass_short.csv")
    tables = [await bacass("bacass_hybrid.csv"), await bacass("bacass_long.csv"), short, short]
    out = await csvtk_concat(tables, in_format="tsv", out_format="csv", prefix="test")
    assert name_of(out) == "test.csv", f"expected test.csv, got {name_of(out)}"
    await assert_md5(out, "bb0ed52999b6b24297bcefb3c29f0a5c", label="test.csv")


tests = [
    test_csvtk_join_csv,
    test_csvtk_join_tsv_output,
    test_csvtk_concat_tsv_to_csv,
]
