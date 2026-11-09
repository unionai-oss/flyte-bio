"""Tests for flyte_bio.modules.gapseq.

The cases are upstream's gapseq tests (modules commit 9a1cc9f, gapseq 2.1.0),
with their chains: find and find-transport feed draft, and draft feeds fill.
As upstream:

- find's tables are checked for their names, line counts and header lines;
- find-transport's table by upstream's "normalised table" md5 (nft-utils'
  ``tableMD5Keys`` with ``tableRemoveComments``: leading comment lines
  dropped, columns sorted by name, decimals rounded to 6 places, absolute
  paths reduced to file names, rows sorted, then written back as CSV), which
  :func:`normed_table_md5` reproduces;
- the models (RDS and SBML, which record dates) by name only;
- the medium by plain md5;
- the database download by its presence.
"""

import csv
import decimal
import hashlib
import io
import math
import os
import re

from flyte.io import File

from flyte_bio.modules.gapseq import (
    gapseq_doall,
    gapseq_draft,
    gapseq_fill,
    gapseq_find,
    gapseq_findtransport,
    gapseq_medium,
    gapseq_requestdb,
)
from tests.framework import assert_md5, assert_nonempty, env, fixture

PORTIERA = "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/"
PROTEOME = PORTIERA + "genome/proteome.fasta"
# The test configs' args.
FIND_ARGS = "-M prot -A diamond -O"
CHAIN_TRANSPORT_ARGS = "-M prot -A diamond"


def name_of(file: File) -> str:
    return os.path.basename(file.path)


def names_of(files: list[File]) -> list[str]:
    return sorted(name_of(f) for f in files)


async def lines_of(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


def java_double(value: float) -> str:
    """``value`` as Java's ``Double.toString`` prints it (shortest round-trip digits)."""
    if value == 0:
        return repr(value)
    if 1e-3 <= abs(value) < 1e7:
        text = f"{decimal.Decimal(repr(value)):f}"
        return text if "." in text else text + ".0"
    digits, exponent = decimal.Decimal(repr(abs(value))).normalize().as_tuple()[1:]
    text = "".join(map(str, digits))
    scientific = len(text) + exponent - 1
    return f"{'-' if value < 0 else ''}{text[0]}.{text[1:] or '0'}E{scientific}"


def is_double(value: str) -> bool:
    if not value.strip() or "_" in value:
        return False
    try:
        float(value)
    except ValueError:
        return False
    return "." in value or "e" in value or "E" in value


def normalize_value(value: str, digits: int) -> str:
    if value.startswith("/") or re.match(r"^[A-Za-z]:[\\/].*", value):
        parts = value.replace("\\", "/").split("/")
        return parts[-1] if len(parts) > 1 else value.replace("\\", "/")
    if is_double(value):
        scale = 10.0**digits
        return java_double(math.floor(float(value) * scale + 0.5) / scale)
    return value


def csv_field(value: str, first: bool) -> str:
    """A field as Apache Commons CSV's DEFAULT format (minimal quoting) prints it."""
    if value == "":
        return '""' if first else ""
    quote = value[0] <= "#" or any(c in value for c in ',"\r\n') or value[-1] <= " "
    if not quote:
        return value
    return '"' + value.replace('"', '""') + '"'


def normed_table_md5(text: str, digits: int = 6, remove_comments: bool = True) -> str:
    """nft-utils' ``md5NormedTable`` of a delimited table."""
    lines = text.splitlines()
    if remove_comments:
        while lines and (not lines[0].strip() or lines[0].strip().startswith("#")):
            lines.pop(0)
        text = "\n".join(lines)
    first = next((line for line in text.splitlines() if line.strip()), "")
    separator = max([",", "\t", ";"], key=lambda s: (first.count(s), -[",", "\t", ";"].index(s)))
    rows = [row for row in csv.reader(io.StringIO(text), delimiter=separator, quotechar='"') if row]
    columns, body = rows[0], rows[1:]
    order = sorted(range(len(columns)), key=lambda i: columns[i])
    columns = [columns[i] for i in order]
    body = [[normalize_value(row[i], digits) for i in order] for row in body]
    body.sort(key="\0".join)
    out = "".join(",".join(csv_field(v, i == 0) for i, v in enumerate(row)) + "\n" for row in [columns, *body])
    return hashlib.md5(out.encode()).hexdigest()


async def find_and_transport(transport_args: str):
    fasta = await fixture(PROTEOME)
    find = await gapseq_find(fasta, args=FIND_ARGS)
    transport = await gapseq_findtransport(fasta, args=transport_args)
    return find, transport


@env.task
async def test_gapseq_requestdb() -> None:
    # upstream case: gapseq/requestdb "bacteria - download database"
    await assert_nonempty(await gapseq_requestdb("Bacteria"), label="gapseq_db.tar.gz")


@env.task
async def test_gapseq_find() -> None:
    # upstream case: gapseq/find "candidatus_portiera_aleyrodidarum proteome [fasta]"
    result = await gapseq_find(await fixture(PROTEOME), args=FIND_ARGS)
    assert name_of(result.reactions) == "proteome-all-Reactions.tbl", name_of(result.reactions)
    assert name_of(result.pathways) == "proteome-all-Pathways.tbl", name_of(result.pathways)
    assert names_of(result.log) == ["proteome-all-find_aligner.log"], names_of(result.log)
    reactions = await lines_of(result.reactions)
    assert len(reactions) > 50, f"reactions: {len(reactions)} lines"
    assert re.match(r"^pathway\trxn\tname\tec", reactions[3]), f"reactions header {reactions[3]!r}"
    pathways = await lines_of(result.pathways)
    assert len(pathways) > 10, f"pathways: {len(pathways)} lines"
    assert re.match(r"^ID\tName", pathways[3]), f"pathways header {pathways[3]!r}"


@env.task
async def test_gapseq_findtransport() -> None:
    # upstream case: gapseq/findtransport "candidatus_portiera_aleyrodidarum proteome [fasta]"
    result = await gapseq_findtransport(await fixture(PROTEOME), args="-M prot -A diamond -b 200")
    assert name_of(result.tbl) == "proteome-Transporter.tbl", name_of(result.tbl)
    assert result.fna == [], names_of(result.fna)
    assert names_of(result.log) == ["proteome-Transporter-find_aligner.log"], names_of(result.log)
    lines = await lines_of(result.tbl)
    assert len(lines) > 5, f"transporters: {len(lines)} lines"
    assert re.match(r"^id\ttc\tsub", lines[3]), f"transporters header {lines[3]!r}"
    actual = normed_table_md5("\n".join(lines) + "\n")
    expected = "a9a30bd83731f7c997dd743aaa85e953"
    assert actual == expected, f"proteome-Transporter.tbl: normalised md5 {actual}, expected {expected}"


@env.task
async def test_gapseq_draft() -> None:
    # upstream case: gapseq/draft "candidatus_portiera_aleyrodidarum proteome [fasta] - real"
    find, transport = await find_and_transport(CHAIN_TRANSPORT_ARGS)
    result = await gapseq_draft(find.reactions, find.pathways, transporters=transport.tbl, prefix="test")
    assert name_of(result.model) == "test-draft.RDS" and name_of(result.xml) == "test-draft.xml"
    assert result.log == [], names_of(result.log)
    await assert_nonempty(result.model, label="test-draft.RDS")
    await assert_nonempty(result.xml, label="test-draft.xml")


@env.task
async def test_gapseq_medium() -> None:
    # upstream case: gapseq/medium "gapseq medium prediction from draft and pathways"
    draft = await fixture(PORTIERA + "gapseq/test-draft.RDS")
    pathways = await fixture(PORTIERA + "gapseq/proteome-all-Pathways.tbl")
    medium = await gapseq_medium(draft, pathways, prefix="test")
    assert name_of(medium) == "test-medium.csv", name_of(medium)
    rows = list(csv.reader(await lines_of(medium)))
    assert len(rows) > 1 and all(len(row) == 3 for row in rows), "medium: expected 3 columns"
    await assert_md5(medium, "c12ef037faf0546ab38d7d195a2f5669", label="test-medium.csv")


@env.task
async def test_gapseq_fill() -> None:
    # upstream case: gapseq/fill "candidatus_portiera_aleyrodidarum proteome [fasta] - real"
    find, transport = await find_and_transport(CHAIN_TRANSPORT_ARGS)
    draft = await gapseq_draft(find.reactions, find.pathways, transporters=transport.tbl, prefix="test")
    result = await gapseq_fill(draft.model, prefix="test")
    assert name_of(result.model) == "test-filled.RDS" and name_of(result.xml) == "test-filled.xml"
    assert result.log == [], names_of(result.log)
    await assert_nonempty(result.model, label="test-filled.RDS")
    await assert_nonempty(result.xml, label="test-filled.xml")


@env.task
async def test_gapseq_doall() -> None:
    # upstream case: gapseq/doall "candidatus_portiera_aleyrodidarum proteome [fasta] + MM_glu medium"
    fasta = await fixture(PROTEOME)
    medium = await fixture("generic/csv/MM_glu.csv")
    result = await gapseq_doall(fasta, medium=medium, args="-A diamond -O")
    assert names_of(result.model) == ["proteome-draft.RDS", "proteome.RDS"], names_of(result.model)
    assert names_of(result.xml) == ["proteome-draft.xml", "proteome.xml"], names_of(result.xml)
    expected_tbl = ["proteome-Transporter.tbl", "proteome-all-Pathways.tbl", "proteome-all-Reactions.tbl"]
    assert names_of(result.tbl) == expected_tbl, names_of(result.tbl)
    expected_log = ["proteome-Transporter-find_aligner.log", "proteome-all-find_aligner.log"]
    assert names_of(result.log) == expected_log, names_of(result.log)
    assert result.fna == [], names_of(result.fna)
    for file in [*result.model, *result.xml, *result.tbl]:
        await assert_nonempty(file, label=name_of(file))


tests = [
    test_gapseq_requestdb,
    test_gapseq_find,
    test_gapseq_findtransport,
    test_gapseq_draft,
    test_gapseq_medium,
    test_gapseq_fill,
    test_gapseq_doall,
]
