"""Tests for flyte_bio.modules.datavzrd.

The cases are upstream's datavzrd tests (modules commit 9a1cc9f, Datavzrd
2.75.3): the report's file list for a CSV, a TSV and a config of both, and
that a run without a config fails.
"""

import asyncio
import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.datavzrd import datavzrd
from tests.framework import env, fixture

CSV = "generic/csv/test.csv"
TSV = "generic/tsv/network.tsv"
NETWORK = [
    "network/config.js",
    "network/data/data_1.js",
    "network/functions.js",
    "network/index_1.html",
    "network/plots/plots.js",
]
TEST = ["test/config.js", "test/data/data_1.js", "test/functions.js", "test/index_1.html", "test/plots/plots.js"]
MULTIPLE_CONFIG = """datasets:
  test:
    path: test.csv
  network:
    path: network.tsv
views:
  test:
    dataset: test
    render-table:
      columns: {}
  network:
    dataset: network
    render-table:
      columns: {}
"""


async def report_files(report: Dir) -> list[str]:
    root = report.path.rstrip("/") + "/"
    return sorted([f.path.removeprefix(root) async for f in report.walk()])


async def check(report: Dir, prefix: str, views: list[str]) -> None:
    assert report.path.rstrip("/").endswith(f"/{prefix}"), report.path
    expected = sorted(["index.html", "static/bundle.js", *views])
    files = await report_files(report)
    assert files == expected, f"report files {files}"


@env.task
async def test_datavzrd_csv() -> None:
    config, csv = await asyncio.gather(fixture("generic/config/config_test.datavzrd.yaml"), fixture(CSV))
    await check(await datavzrd(config, [csv], prefix="csv_config"), "csv_config", TEST)


@env.task
async def test_datavzrd_tsv() -> None:
    config, tsv = await asyncio.gather(fixture("generic/config/config_network.datavzrd.yaml"), fixture(TSV))
    await check(await datavzrd(config, [tsv], prefix="tsv_config"), "tsv_config", NETWORK)


@env.task
async def test_datavzrd_multiple_tables() -> None:
    path = Path(tempfile.mkdtemp()) / "config_multiple.datavzrd.yaml"
    path.write_text(MULTIPLE_CONFIG)
    config = await File.from_local(str(path))
    files = list(await asyncio.gather(fixture(CSV), fixture(TSV)))
    await check(await datavzrd(config, files, prefix="multiple_config"), "multiple_config", NETWORK + TEST)


@env.task
async def test_datavzrd_no_config_file() -> None:
    # As upstream, the run fails; its failure arrives as a generic pod error.
    try:
        await datavzrd(None, [], prefix="test")
    except Exception:
        return
    raise AssertionError("datavzrd ran without a config file")


tests = [test_datavzrd_csv, test_datavzrd_tsv, test_datavzrd_multiple_tables, test_datavzrd_no_config_file]
