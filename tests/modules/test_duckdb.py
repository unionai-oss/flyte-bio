"""Tests for flyte_bio.modules.duckdb.

The cases and expected md5s are upstream's duckdb/table2parquet tests (modules
commit 9a1cc9f, DuckDB 1.5.5), with the ``ext.args`` and ``ext.prefix`` their
config sets: ``nullstr=['NA','']`` for the ``nullstr-*`` tables, and the
heredoc-forging prefix. That prefix holds newlines, which Flyte's object store
can't hold in a file name, so its case checks that the task runs it and that
the wrapper refuses it rather than the output's name. The gzipped cases compress the fixture with pigz first,
as upstream. The generated tables are built here as upstream's test builds them.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.duckdb import duckdb_table2parquet, duckdb_table2parquet_cmd
from flyte_bio.modules.pigz import pigz_compress
from tests.framework import assert_md5, env, fixture

CSV = "generic/csv/test.csv"
TSV = "generic/tsv/test.tsv"
NULLSTR = "nullstr=['NA','']"


def name_of(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def local_table(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


def na_table(rows: int, is_na) -> str:
    """Upstream's generated two-column table: ``g<i>`` and ``i % 7``, or ``NA`` where ``is_na(i)``."""
    return "id\tvalue\n" + "".join(f"g{i}\t{'NA' if is_na(i) else i % 7}\n" for i in range(1, rows + 1))


async def parquet_contains_na(file: File) -> bool:
    async with file.open("rb") as fh:
        return b"NA" in bytes(await fh.read())


@env.task
async def test_duckdb_csv() -> None:
    out = await duckdb_table2parquet(await fixture(CSV), prefix="test")
    assert name_of(out) == "test.parquet", f"expected test.parquet, got {name_of(out)}"
    await assert_md5(out, "62f39a0ec6f40def92860479260b4e16", label="test.parquet")


@env.task
async def test_duckdb_csv_gzipped() -> None:
    out = await duckdb_table2parquet(await pigz_compress(await fixture(CSV)), prefix="test")
    await assert_md5(out, "62f39a0ec6f40def92860479260b4e16", label="test.parquet")


@env.task
async def test_duckdb_tsv() -> None:
    out = await duckdb_table2parquet(await fixture(TSV), prefix="test")
    await assert_md5(out, "180fa2a3c593e6486a26180242211b27", label="test.parquet")


@env.task
async def test_duckdb_tsv_gzipped() -> None:
    out = await duckdb_table2parquet(await pigz_compress(await fixture(TSV)), prefix="test")
    await assert_md5(out, "180fa2a3c593e6486a26180242211b27", label="test.parquet")


@env.task
async def test_duckdb_prefix_with_quote_and_shell_characters() -> None:
    prefix = "O'Brien$(touch pwned)`echo also_pwned`"
    out = await duckdb_table2parquet(await fixture(CSV), prefix=prefix)
    assert name_of(out) == f"{prefix}.parquet", f"expected {prefix}.parquet, got {name_of(out)}"


@env.task
async def test_duckdb_prefix_with_newline_and_heredoc_delimiter() -> None:
    # Flyte's object store can't address a name containing a newline, so the
    # wrapper refuses such a prefix; the task itself still runs it safely.
    prefix = "safe\nSQL\ntouch pwned"
    table = await fixture(CSV)
    await duckdb_table2parquet_cmd(table=table, table_name="csv", delimiter=",", prefix=prefix)
    try:
        await duckdb_table2parquet(table, prefix=prefix)
    except ValueError as e:
        assert "newline" in str(e), str(e)
    else:
        raise AssertionError("duckdb_table2parquet accepted a prefix with a newline")


@env.task
async def test_duckdb_nullstr_after_type_detection_sample() -> None:
    # DuckDB types each column from the first 20480 rows; a numeric column holds no "NA" text.
    table = await local_table("late_na.tsv", na_table(25000, lambda i: i > 24000))
    out = await duckdb_table2parquet(table, prefix="nullstr-late", read_options=NULLSTR)
    assert not await parquet_contains_na(out), "the value column was typed as text"
    await assert_md5(out, "817002b8a6daece09529c9f1671de210", label="nullstr-late.parquet")


@env.task
async def test_duckdb_nullstr_in_first_rows() -> None:
    table = await local_table("early_na.tsv", na_table(30, lambda i: i <= 3))
    out = await duckdb_table2parquet(table, prefix="nullstr-early", read_options=NULLSTR)
    assert not await parquet_contains_na(out), "the value column was typed as text"
    await assert_md5(out, "7ababb9d8a62166773b0bc071cc5ec35", label="nullstr-early.parquet")


@env.task
async def test_duckdb_na_is_text_without_nullstr() -> None:
    table = await local_table("text_na.tsv", na_table(30, lambda i: i <= 3))
    out = await duckdb_table2parquet(table, prefix="plain-na")
    assert await parquet_contains_na(out), "expected NA to be kept as text"
    await assert_md5(out, "dc828d68b9b16698c5f5cc8fc191d42c", label="plain-na.parquet")


tests = [
    test_duckdb_csv,
    test_duckdb_csv_gzipped,
    test_duckdb_tsv,
    test_duckdb_tsv_gzipped,
    test_duckdb_prefix_with_quote_and_shell_characters,
    test_duckdb_prefix_with_newline_and_heredoc_delimiter,
    test_duckdb_nullstr_after_type_detection_sample,
    test_duckdb_nullstr_in_first_rows,
    test_duckdb_na_is_text_without_nullstr,
]
