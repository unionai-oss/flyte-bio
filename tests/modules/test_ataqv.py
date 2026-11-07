"""Tests for flyte_bio.modules.ataqv.

The cases and expected values are upstream's ataqv/ataqv and ataqv/mkarv
tests (modules commit 9a1cc9f, ataqv 1.3.1), all with organism ``human``.

ataqv/ataqv snapshots the md5 of Groovy's ``toString()`` of the metrics
file's first record without its ``timestamp``; :func:`groovy_string`
renders a parsed record the same way (maps as ``[key:value, ...]``, lists as
``[a, b]``, JSON decimals as Java's ``BigDecimal`` prints them, which
Python's ``Decimal`` matches). ataqv/mkarv snapshots every report file's md5,
the data file's decompressed.
"""

import asyncio
import decimal
import hashlib
import json

from flyte.io import File

from flyte_bio.modules.ataqv import ataqv_ataqv, ataqv_mkarv
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

BAM = "genomics/sarscov2/illumina/bam/"
BED = "genomics/sarscov2/genome/bed/"
# The record's md5 without and with a TSS file; the peak and excluded-region files don't change it.
PLAIN_MD5 = "4d418b039281d274b1440c9e73267143"
TSS_MD5 = "4fb5f3438d800d24f5019d232aef082d"

REPORT_MD5S = {
    "css/ataqv.css": "7bb9e18974018dd29bdc12a46a6d6d65",
    "css/datatables.buttons.min.css": "675306abc9d3239314f2c66c0e22ef7a",
    "css/datatables.fontawesome.css": "716d7a3edd2c75e712b1ffe26f12a861",
    "css/datatables.min.css": "00ef3e39d9302618395ab189ec8bcd9f",
    "css/font-awesome.min.css": "89916fa773ce96569604016ef25cab50",
    "css/normalize.css": "40d20cc85100840580df4346bb99d7ed",
    "fonts/FontAwesome.otf": "9ccfa32dd4cd1b8e83f68899d85bd5e6",
    "fonts/fontawesome-webfont.eot": "404a525502f8e5ba7e93b9f02d9e83a9",
    "fonts/fontawesome-webfont.svg": "bae4a87c1e5dff40baa3f49d52f5347a",
    "fonts/fontawesome-webfont.ttf": "fb650aaf10736ffb9c4173079616bf01",
    "fonts/fontawesome-webfont.woff": "891e3f340c1126b4c7c142e5f6e86816",
    "fonts/fontawesome-webfont.woff2": "926c93d201fe51c8f351e858468980c3",
    "fonts/sourcesanspro-regular.woff": "f7bd788f18b8c4bb93dd37d140348e1e",
    "fonts/sourcesanspro-regularit.woff": "c3638b17f4fd76dae12fe2ae14571e57",
    "fonts/sourcesanspro-semibold.woff": "e7fc8925d9364e9d177d9e1d08bb1855",
    "fonts/sourcesanspro-semiboldit.woff": "f1d255aa459786dfc6aa2e488ac01245",
    "index.html": "bf7747be761e56ad7c54c842ac88461a",
    "js/ataqv.js": "feb291b7839e9e43ed304565e3a605d9",
    "js/configuration.js": "073bc07f1371e28a97098f8b10469ae4",
    "js/d3.min.js": "db69fb2626a71a286ee772d673138aca",
    "js/datatables.min.js": "e369b872620dadb05e4eb555b81f9112",
    "js/jszip.min.js": "09e492cb492ffa75484bbe10f1f721d1",
}
DATA_MD5 = "7c01e4f0136404e6db2328485d68d154"  # data/1.json.gz, decompressed


def groovy_string(value: object) -> str:
    """``value`` as Groovy's ``toString()`` prints the map JsonSlurper parses it to."""
    if isinstance(value, dict):
        if not value:
            return "[:]"
        return "[" + ", ".join(f"{key}:{groovy_string(item)}" for key, item in value.items()) + "]"
    if isinstance(value, list):
        return "[" + ", ".join(groovy_string(item) for item in value) + "]"
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


async def record_md5(metrics: File) -> str:
    async with metrics.open("rb") as fh:
        records = json.loads(bytes(await fh.read()), parse_float=decimal.Decimal)
    record = {key: value for key, value in records[0].items() if key != "timestamp"}
    return hashlib.md5(groovy_string(record).encode()).hexdigest()


async def check_metrics(metrics: File, expected: str, label: str) -> None:
    assert metrics.name == "test.ataqv.json", f"{label}: metrics file is {metrics.name}"
    actual = await record_md5(metrics)
    assert actual == expected, f"{label}: record md5 mismatch (expected {expected}, got {actual})"


@env.task
async def test_ataqv_ataqv() -> None:
    # upstream cases: ataqv/ataqv "test_ataqv_ataqv" and "test_ataqv_ataqv_problem_reads"
    bam = await fixture(BAM + "test.paired_end.bam")
    plain, problems = await asyncio.gather(
        ataqv_ataqv(bam, "human", prefix="test"),
        ataqv_ataqv(bam, "human", prefix="test", args="--log-problematic-reads"),
    )
    await check_metrics(plain.json, PLAIN_MD5, "plain")
    assert plain.problems == [], f"plain: wrote {[f.name for f in plain.problems]}"
    await check_metrics(problems.json, PLAIN_MD5, "problem reads")
    assert [f.name for f in problems.problems] == ["1.problems"], [f.name for f in problems.problems]


@env.task
async def test_ataqv_ataqv_regions() -> None:
    # upstream cases: ataqv/ataqv "test_ataqv_ataqv_peak", "test_ataqv_ataqv_tss" and "test_ataqv_ataqv_excluded_regs"
    bam, bai, bed, bed2 = await asyncio.gather(
        fixture(BAM + "test.paired_end.sorted.bam"),
        fixture(BAM + "test.paired_end.sorted.bam.bai"),
        fixture(BED + "test.bed"),
        fixture(BED + "test2.bed"),
    )
    peak, tss, excluded = await asyncio.gather(
        ataqv_ataqv(bam, "human", prefix="test", peak_file=bed),
        ataqv_ataqv(bam, "human", prefix="test", bai=bai, tss_file=bed),
        ataqv_ataqv(bam, "human", prefix="test", bai=bai, tss_file=bed, excluded_regions_file=bed2),
    )
    for label, result, md5 in (("peak", peak, PLAIN_MD5), ("tss", tss, TSS_MD5), ("excluded", excluded, TSS_MD5)):
        await check_metrics(result.json, md5, label)
        assert result.problems == [], f"{label}: wrote {[f.name for f in result.problems]}"


@env.task
async def test_ataqv_mkarv() -> None:
    # upstream case: ataqv/mkarv "test_ataqv_mkarv"
    metrics = await ataqv_ataqv(await fixture(BAM + "test.paired_end.bam"), "human", prefix="test")
    html = await ataqv_mkarv([metrics.json])
    base = html.path.rstrip("/")
    files = {file.path[len(base) :].lstrip("/"): file async for file in html.walk()}
    assert sorted(files) == sorted([*REPORT_MD5S, "data/1.json.gz"]), f"report files: {sorted(files)}"
    for name, md5 in REPORT_MD5S.items():
        await assert_md5(files[name], md5, label=name)
    await assert_gunzipped_md5(files["data/1.json.gz"], DATA_MD5, label="data/1.json.gz")


tests = [test_ataqv_ataqv, test_ataqv_ataqv_regions, test_ataqv_mkarv]
