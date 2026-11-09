"""Tests for flyte_bio.modules.hmmer.

The cases are upstream's hmmer/hmmsearch, hmmer/hmmalign, hmmer/hmmbuild and
hmmer/hmmrank tests (modules commit 9a1cc9f, HMMER 3.4, DuckDB 1.5.5).

- hmmsearch: the bacterial 16S rRNA profile searched against E. coli K-12's
  16S gene. Upstream snapshots only the alignment's md5 (decompressed, as
  nf-test hashes ``.gz`` outputs); the report ends with run-time statistics and
  the tables record the run's directory and date, so those get content checks.
- hmmalign: the same gene aligned to the same profile, by md5.
- hmmbuild: upstream marks both outputs unstable (they record the run's date
  and CPU time), so they get content checks.
- hmmrank: upstream's setup fetches the 16S profiles from barrnap's archaeal
  and bacterial databases, searches them against the phyloplace 16S test
  sequences, converts the tables to TSV with its HMMER table formatter, and
  then to Parquet. The fetch and the formatter aren't modules here, so
  :func:`fetch_profile` and :func:`format_table` do the same in Python (the
  formatter's parsing rules, unchanged); the search and the Parquet conversion
  use :func:`hmmsearch` and :func:`duckdb_table2parquet`. The curated NrdJ
  tables and the inline NA and tied-start tables are upstream's. The
  heredoc-forging prefix holds newlines, which Flyte's object store can't hold
  in a file name, so its case checks that the task runs it and that the
  wrapper refuses it, rather than the output's name.
- esl_sfetch_index: upstream's hmmer/eslsfetchindex test (Easel 0.49), the
  SARS-CoV-2 scaffolds' SSI index, by md5.
- esl_alimask and esl_reformat: upstream's hmmer/eslalimask and
  hmmer/eslreformat tests (Easel 0.49), on hmmalign's alignment of the 16S
  gene (upstream's setup), by md5 (decompressed). Upstream's unalign case pipes
  the output through ``sed "/^>/!s/-//g"``; here that's the ``sed_script``.
- format_tsv: upstream's hmmer/formattsv tests (gawk 5.3.1). Where upstream
  snapshots the TSV's lines, the TSV is checked by the md5 of those lines
  (each newline-terminated, as the formatter writes them), with upstream's
  row-count assertions. The synthetic tables are upstream's, written under
  their temp files' prefixes. In the CRLF case upstream's table ends in
  ``\\r\\r\\n``, so the formatter leaves one carriage return on the last
  row, which upstream's assertions (Java's ``readLines``, which ends a line at
  a carriage return too) don't see; the lines are split the same way here. The
  failing cases check that the task fails; when the error carries the
  formatter's message, that it's upstream's.
"""

import asyncio
import gzip
import re
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.duckdb import duckdb_table2parquet
from flyte_bio.modules.hmmer import (
    esl_alimask,
    esl_reformat,
    esl_sfetch_index,
    format_tsv,
    hmmalign,
    hmmbuild,
    hmmrank,
    hmmrank_cmd,
    hmmsearch,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HMM = "delete_me/hmmer/bac.16S_rRNA.hmm.gz"
SEQDB = "delete_me/hmmer/e_coli_k12_16s.fna.gz"
TARGET = "CP025268.1:223601-225154"


async def gunzipped_text(file: File) -> str:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode()


async def check_report(report: File) -> None:
    text = await gunzipped_text(report)
    assert "Query:       16S_rRNA" in text, "hmmsearch report lacks the query"
    assert TARGET in text, "hmmsearch report lacks the 16S hit"
    assert text.rstrip().endswith("[ok]"), "hmmsearch report is incomplete"


async def check_table(table: File, label: str) -> None:
    rows = [line.split() for line in (await gunzipped_text(table)).splitlines() if not line.startswith("#")]
    assert rows, f"{label} has no hits"
    assert all(row[0] == TARGET for row in rows), f"{label} hits {sorted({row[0] for row in rows})}"
    assert all("16S_rRNA" in row for row in rows), f"{label} rows lack the query name"


@env.task
async def test_hmmsearch_txt() -> None:
    result = await hmmsearch(await fixture(HMM), await fixture(SEQDB), prefix="test")
    assert result.output.path.endswith("/test.txt.gz")
    await check_report(result.output)
    assert result.alignments is None and result.target_summary is None and result.domain_summary is None
    assert sorted(f.name for f in await result.results.list_files()) == ["test.txt.gz"]


@env.task
async def test_hmmsearch_txt_sto_tbl_domtbl() -> None:
    result = await hmmsearch(
        await fixture(HMM), await fixture(SEQDB), prefix="test", write_align=True, write_target=True, write_domain=True
    )
    await check_report(result.output)
    assert result.alignments is not None and result.target_summary is not None and result.domain_summary is not None
    await assert_gunzipped_md5(result.alignments, "5c44c289b9e36aa1f7f3afae2005fbb7", label="test.sto.gz")
    await check_table(result.target_summary, "test.tbl.gz")
    await check_table(result.domain_summary, "test.domtbl.gz")


@env.task
async def test_hmmalign() -> None:
    out = await hmmalign(await fixture(SEQDB), await fixture(HMM), prefix="test")
    assert out.path.endswith("/test.sto.gz"), out.path
    await assert_gunzipped_md5(out, "4ae989d5ade2aaae9578cb88ba031e8f", label="test.sto.gz")


@env.task
async def test_hmmbuild() -> None:
    result = await hmmbuild(await fixture("delete_me/hmmer/PF14720_seed.alnfaa.gz"), prefix="PF14720")
    assert result.hmm.path.endswith("/PF14720.hmm.gz"), result.hmm.path
    assert result.report.path.endswith("/PF14720.hmmbuild.txt"), result.report.path
    profile = await gunzipped_text(result.hmm)
    assert profile.startswith("HMMER3/f") and "\nNAME  PF14720\n" in profile, "hmmbuild wrote no PF14720 profile"
    async with result.report.open("rb") as fh:
        report = bytes(await fh.read()).decode()
    assert "CPU time:" in report, "hmmbuild report lacks its CPU time"


# hmmrank's setup: upstream's inputs, and the formatter's TSV columns.
PHYLOPLACE = "https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/"
BARRNAP = "https://raw.githubusercontent.com/tseemann/barrnap/0.9/db/"
NRDJ_TABLES = {"tblout": ["NrdJ.tbl.gz", "NrdJm1.tbl.gz"], "domtblout": ["NrdJ.domtbl.gz", "NrdJm1.domtbl.gz"]}
TABLE_COLUMNS = {
    "tblout": "profile target_name target_accession query_name query_accession full_evalue full_score full_bias "
    "best_domain_evalue best_domain_score best_domain_bias exp reg clu ov env dom rep inc description",
    "domtblout": "profile target_name target_accession target_length query_name query_accession query_length "
    "full_evalue full_score full_bias domain_number domain_count c_evalue i_evalue domain_score domain_bias "
    "hmm_from hmm_to ali_from ali_to env_from env_to acc description",
}
COORDINATE_COLUMNS = [
    "tlen",
    "qlen",
    "hmm_from",
    "hmm_to",
    "hmm_len",
    "hmm_n_islands",
    "ali_from",
    "ali_to",
    "ali_len",
    "ali_n_islands",
    "env_from",
    "env_to",
    "env_len",
    "env_n_islands",
]

# Upstream's inline tables: seq2 is a hit with no domain rows; seq1's tied-start domains merge into two islands.
NA_TBLOUT = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 has a domain match
seq2                               -          model_A               PF00001       2e-10   49.0   0.1      3e-10   48.0   0.1   1.0   1   0   0   1   1   1   1 no domain match
"""  # noqa: E501
NA_DOMTBLOUT = """\
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen  E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
seq1                               -   100 model_A               PF00001   100      1e-10   50.0   0.1   1   1     1e-10      1e-10   50.0   0.1     1    50     1    50     1    50 0.95 has a domain match
"""  # noqa: E501
TIE_TBLOUT = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   3.0   3   0   0   3   3   3   1 tied domain starts
"""  # noqa: E501
TIE_DOMTBLOUT = """\
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen  E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
seq1                               -   200 model_A               PF00001   200      1e-10   50.0   0.1   1   3     1e-10      1e-10   50.0   0.1     1     3     1     3     1     3 0.95 first domain
seq1                               -   200 model_A               PF00001   200      1e-10   50.0   0.1   2   3     1e-10      1e-10   50.0   0.1     5   100     5   100     5   100 0.95 tied start, far end
seq1                               -   200 model_A               PF00001   200      1e-10   50.0   0.1   3   3     1e-10      1e-10   50.0   0.1     5     6     5     6     5     6 0.95 tied start, near end
"""  # noqa: E501

# One row per (target, profile or model): accno, profile, the hmm/ali/env bounds, lengths and island
# counts, then tlen and qlen.
NRDJ_EXPECTED = sorted(
    " ".join(row.split())
    for row in [
        "AACBMDAH_04658 NrdJ   11  95  68 2  57  238  67 2  40  242 128 2  805 407",
        "AACBMDAH_04658 NrdJm1  1 397 397 1  14  803 400 2  14  805 425 2  805 748",
        "AACBMDAH_04659 NrdJm1 513 748 236 1   1  620 239 2   1  620 249 2  626 748",
        "AAGAFPGM_00302 NrdJ   12 406 379 3 176  724 377 5 156  725 489 5  780 407",
        "AAGAFPGM_00302 NrdJm1  6 711 588 4 134  735 555 4 129  770 615 2  780 748",
        "DPHPGHLN_02900 NrdJ   62 384 306 2 205 1379 327 4 190 1397 401 4 1460 407",
        "DPHPGHLN_02900 NrdJm1  1 747 747 1  14 1454 752 3  14 1455 791 3 1460 748",
        "MNMMFBCN_00658 NrdJ   11 386 353 3  84  718 374 5  67  732 485 4  798 407",
        "MNMMFBCN_00658 NrdJm1  1 748 748 1  42  792 751 1  42  792 751 1  798 748",
    ]
)
NRDJ_COLUMNS = COORDINATE_COLUMNS[2:] + COORDINATE_COLUMNS[:2]


def is_dashes(token: str) -> bool:
    return re.fullmatch(r"-+", token) is not None


def format_rows(text: str, label: str, table_format: str) -> list[str]:
    """One HMMER table's rows as TSV lines, labelled ``label``: upstream's formatter, in Python.

    The header row of dashes gives the column count N; each data row's first
    N-1 fields are columns and the rest is rejoined into the description. A
    data row after the run-metadata footer (a bare ``#``) is an error.
    """
    expected = len(TABLE_COLUMNS[table_format].split()) - 1
    lines = text.split("\n")
    if text.endswith("\n"):
        lines.pop()
    columns, in_footer, rows = 0, False, []
    for raw in lines:
        line = raw.removesuffix("\r")
        fields = [f for f in re.split(r"[ \t]+", line) if f]
        if columns == 0 and line.startswith("#"):
            first = fields[0][1:]
            count, ok = 0, True
            if first:
                ok = is_dashes(first)
                count += ok
            if ok:
                for field in fields[1:]:
                    if not is_dashes(field):
                        ok = False
                        break
                    count += 1
            if ok and count > 0:
                if count != expected:
                    raise ValueError(f"expected {expected} {table_format} columns, found {count}")
                columns = count
        elif columns > 0 and line == "#":
            in_footer = True
        elif in_footer:
            if line and not line.startswith("#"):
                raise ValueError("data row seen after the run-metadata footer")
        elif columns == 0:
            raise ValueError("data row seen before a column-count header was found")
        else:
            padded = fields + [""] * (columns - len(fields))
            rows.append("\t".join([label, *padded[: columns - 1], " ".join([padded[columns - 1], *fields[columns:]])]))
    return rows


def unzipped(data: bytes) -> str:
    return (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()


async def file_bytes(file: File) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read())


def local_file(name: str, data: bytes) -> str:
    path = Path(tempfile.mkdtemp()) / name
    path.write_bytes(data)
    return str(path)


async def format_table(case: str, labels: list[str], texts: list[str], table_format: str) -> File:
    """``<case>.<format>.tsv.gz``: the tables' rows under one header, as upstream's formatter writes it."""
    rows = [row for label, text in zip(labels, texts) for row in format_rows(text, label, table_format)]
    tsv = "\t".join(TABLE_COLUMNS[table_format].split()) + "\n" + "".join(row + "\n" for row in rows)
    return await File.from_local(local_file(f"{case}.{table_format}.tsv.gz", gzip.compress(tsv.encode(), mtime=0)))


def fetch_profile(hmm_text: str, name: str) -> str:
    """The profile called ``name`` from a multi-profile HMM file, as ``hmmfetch <file> <name>`` selects it."""
    for record in hmm_text.split("//\n"):
        if re.search(rf"^NAME +{re.escape(name)}$", record, re.MULTILINE):
            return record.lstrip("\n") + "//\n"
    raise ValueError(f"no profile {name}")


def multimodel_table(texts: list[str]) -> str:
    """Upstream's several-models table: the first table's comment lines, then every table's data rows."""
    header = [line for line in texts[0].splitlines() if line.startswith("#")]
    data = [line for text in texts for line in text.splitlines() if not line.startswith("#")]
    return "\n".join(header + data)


async def search_16s(database: str) -> tuple[str, str]:
    """The 16S profile of a barrnap database searched against the 16S test sequences: (tblout, domtblout)."""
    profile = fetch_profile(unzipped(await file_bytes(await fixture(BARRNAP + f"{database}.hmm"))), "16S_rRNA")
    hmm = await File.from_local(local_file(f"{database}16s.hmm", profile.encode()))
    result = await hmmsearch(
        hmm, await fixture(PHYLOPLACE + "domain_16s.fna"), prefix=f"{database}16s", write_target=True, write_domain=True
    )
    assert result.target_summary is not None and result.domain_summary is not None
    return unzipped(await file_bytes(result.target_summary)), unzipped(await file_bytes(result.domain_summary))


@env.task(cache="auto")
async def hmmrank_parquet_tables(case: str) -> list[File]:
    """Upstream's setup for one hmmrank case: its tblout and domtblout tables, as Parquet."""
    if case == "16S-test":
        searches = await asyncio.gather(search_16s("arc"), search_16s("bac"))
        labels = ["arc16s", "bac16s"]
        texts = {"tblout": [s[0] for s in searches], "domtblout": [s[1] for s in searches]}
    elif case in ("nrdj", "nrdjboth"):
        texts = {
            fmt: [unzipped(await file_bytes(await fixture(PHYLOPLACE + "hmmrank/" + n))) for n in names]
            for fmt, names in NRDJ_TABLES.items()
        }
        labels = ["NrdJ", "NrdJm1"]
        if case == "nrdjboth":
            texts = {fmt: [multimodel_table(t)] for fmt, t in texts.items()}
            labels = ["nrdj"]
    elif case == "na-test":
        labels, texts = ["P"], {"tblout": [NA_TBLOUT], "domtblout": [NA_DOMTBLOUT]}
    elif case == "tie-test":
        labels, texts = ["P"], {"tblout": [TIE_TBLOUT], "domtblout": [TIE_DOMTBLOUT]}
    else:
        raise ValueError(f"no hmmrank case {case}")
    return [
        await duckdb_table2parquet(await format_table(case, labels, texts[fmt], fmt), prefix=f"{case}.{fmt}")
        for fmt in ("tblout", "domtblout")
    ]


async def ranked_rows(ranks: File) -> tuple[list[str], list[list[str]]]:
    lines = unzipped(await file_bytes(ranks)).splitlines()
    return lines[0].split("\t"), [line.split("\t") for line in lines[1:]]


async def check_known_table(ranks: File, key: str) -> None:
    header, rows = await ranked_rows(ranks)
    columns = ["accno", key, *NRDJ_COLUMNS]
    got = sorted(" ".join(row[header.index(c)] for c in columns) for row in rows)
    assert got == NRDJ_EXPECTED, f"hmmrank rows {got}"


@env.task
async def test_hmmrank_16s() -> None:
    tblout, _ = await hmmrank_parquet_tables("16S-test")
    ranks = await hmmrank(tblout, prefix="16S-test")
    assert ranks.path.endswith("/16S-test.hmmrank.tsv.gz"), ranks.path
    await assert_gunzipped_md5(ranks, "76f4f8aa2d63cb35e81ab31e72260f67", label="16S-test.hmmrank.tsv.gz")


@env.task
async def test_hmmrank_16s_domain_coordinates() -> None:
    tblout, domtblout = await hmmrank_parquet_tables("16S-test")
    ranks = await hmmrank(tblout, domtblout, prefix="16S-test")
    header, rows = await ranked_rows(ranks)
    assert set(COORDINATE_COLUMNS) <= set(header), f"hmmrank header {header}"
    assert rows, "hmmrank wrote no rows"
    # NA is legitimate (a hit with no domain rows), but all-NA is what a broken join looks like.
    lengths = [row[header.index("hmm_len")] for row in rows]
    assert all(v.isdigit() or v == "NA" for v in lengths) and any(v.isdigit() for v in lengths), lengths
    await assert_gunzipped_md5(ranks, "d0cf9473c7a84ae0b65345ad0e1b3fd1", label="16S-test.hmmrank.tsv.gz")


@env.task
async def test_hmmrank_known_domain_coordinates() -> None:
    tblout, domtblout = await hmmrank_parquet_tables("nrdj")
    await check_known_table(await hmmrank(tblout, domtblout, prefix="nrdj"), "profile")


@env.task
async def test_hmmrank_several_models_in_one_table() -> None:
    tblout, domtblout = await hmmrank_parquet_tables("nrdjboth")
    ranks = await hmmrank(tblout, domtblout, prefix="nrdjboth")
    _, rows = await ranked_rows(ranks)
    assert len(rows) == len(NRDJ_EXPECTED), f"expected one row per (target, model), got {len(rows)}"
    await check_known_table(ranks, "profile_desc")


@env.task
async def test_hmmrank_hit_without_domains_is_na() -> None:
    tblout, domtblout = await hmmrank_parquet_tables("na-test")
    header, rows = await ranked_rows(await hmmrank(tblout, domtblout, prefix="na-test"))
    by_accno = {row[header.index("accno")]: row for row in rows}
    assert by_accno["seq1"][header.index("hmm_len")] == "50", by_accno["seq1"]
    assert by_accno["seq2"][header.index("hmm_len")] == "NA", by_accno["seq2"]


@env.task
async def test_hmmrank_tied_domain_starts() -> None:
    tblout, domtblout = await hmmrank_parquet_tables("tie-test")
    header, rows = await ranked_rows(await hmmrank(tblout, domtblout, prefix="tie-test"))
    got = {c: rows[0][header.index(c)] for c in ("hmm_from", "hmm_to", "hmm_len", "hmm_n_islands")}
    assert got == {"hmm_from": "1", "hmm_to": "100", "hmm_len": "99", "hmm_n_islands": "2"}, got


async def check_prefix(prefix: str) -> None:
    tblout, _ = await hmmrank_parquet_tables("nrdj")
    ranks = await hmmrank(tblout, prefix=prefix)
    name = ranks.path.rstrip("/").rsplit("/", 1)[-1]
    assert name == f"{prefix}.hmmrank.tsv.gz", f"expected {prefix!r}.hmmrank.tsv.gz, got {name!r}"
    _, rows = await ranked_rows(ranks)
    assert rows, "hmmrank wrote no rows"


@env.task
async def test_hmmrank_prefix_with_quote() -> None:
    await check_prefix("O'Brien")


@env.task
async def test_hmmrank_prefix_with_shell_characters() -> None:
    await check_prefix("test$(touch pwned)`echo also_pwned`")


@env.task
async def test_hmmrank_prefix_with_newline_and_heredoc_delimiter() -> None:
    # Flyte's object store can't address a name containing a newline, so the
    # wrapper refuses such a prefix; the task itself still runs it safely.
    prefix = "safe\nSQL\ntouch pwned"
    tblout, _ = await hmmrank_parquet_tables("nrdj")
    await hmmrank_cmd(tblout=tblout, domtblout=None, prefix=prefix)
    try:
        await hmmrank(tblout, prefix=prefix)
    except ValueError as e:
        assert "newline" in str(e), str(e)
    else:
        raise AssertionError("hmmrank accepted a prefix with a newline")


@env.task
async def test_hmmer_esl_sfetch_index() -> None:
    ssi = await esl_sfetch_index(await fixture("genomics/sarscov2/illumina/fasta/scaffolds.fasta"))
    assert ssi.name == "scaffolds.fasta.ssi", ssi.name
    await assert_md5(ssi, "759746a03598daa237d86bd6f5030859", label="scaffolds.fasta.ssi")


async def aligned_16s() -> File:
    """Upstream's setup for esl-alimask and esl-reformat: the 16S gene aligned to the 16S profile."""
    return await hmmalign(await fixture(SEQDB), await fixture(HMM), prefix="test")


@env.task
async def test_hmmer_esl_alimask_no_masks() -> None:
    result = await esl_alimask(await aligned_16s(), prefix="test", args="--rf-is-mask")
    assert result.masked.name == "test.masked.sthlm.gz", result.masked.name
    assert result.masks == {}, sorted(result.masks)
    await assert_gunzipped_md5(result.masked, "9e2d2c1f4dbe2851457b478fd75e7f69", label="test.masked.sthlm.gz")


ALIMASK_MD5S = {
    "fmask-rf": "9ec8400b91b0ad4b69f91888f6762caf",
    "fmask-all": "7da678b3d106269ef4675b575e313244",
    "gmask-rf": "5e028068b40861990b5e8c18cb4fbb75",
    "gmask-all": "e1a8c3410bd77bb4c27410c8ec5f35ea",
    "pmask-rf": "9ec8400b91b0ad4b69f91888f6762caf",
    "pmask-all": "7da678b3d106269ef4675b575e313244",
}


@env.task
async def test_hmmer_esl_alimask_all_masks() -> None:
    flags = {f"write_{m.replace('-', '_')}": True for m in ALIMASK_MD5S}
    result = await esl_alimask(await aligned_16s(), prefix="test", args="-pg", **flags)
    await assert_gunzipped_md5(result.masked, "5313874ae46c8a913f507c507321a3d9", label="test.masked.sthlm.gz")
    assert sorted(result.masks) == sorted(ALIMASK_MD5S), sorted(result.masks)
    for mask, md5 in ALIMASK_MD5S.items():
        assert result.masks[mask].name == f"test.{mask}.gz", result.masks[mask].name
        await assert_gunzipped_md5(result.masks[mask], md5, label=f"test.{mask}.gz")


async def check_reformat(args: str, name: str, md5: str, sed_script: str = "") -> None:
    out = await esl_reformat(await aligned_16s(), prefix="test", sed_script=sed_script, args=args)
    assert out.name == name, out.name
    await assert_gunzipped_md5(out, md5, label=name)


@env.task
async def test_hmmer_esl_reformat_afa() -> None:
    await check_reformat("afa", "test.afa.gz", "602548f6a1633eb17a674dfeb43f4724")


@env.task
async def test_hmmer_esl_reformat_phylip() -> None:
    await check_reformat("-u phylip", "test.phylip.gz", "c5122abb9e54b4b6dc8090c05eaaabed")


@env.task
async def test_hmmer_esl_reformat_unalign() -> None:
    await check_reformat("--gapsym=- afa", "test.afa.gz", "f6e80a6bb1ed8eef4cfba4e1980e649a", sed_script="/^>/!s/-//g")


# Upstream's hmmer/formattsv tables.
TSV_TBLOUT_A = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 a multi word description
seq2                               -          model_A               PF00001       3e-08   40.2   0.0      4e-08   39.9   0.0   1.0   1   0   0   1   1   1   1 another one
"""  # noqa: E501
TSV_TBLOUT_B = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                       accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ---------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq3_with_a_longer_name              -          model_B               PF00002       5e-20   99.9   1.2      6e-20   99.0   1.2   1.0   1   0   0   1   1   1   1 third hit, single word desc
"""  # noqa: E501
TSV_DOMTBLOUT_C = """\
#                                                                                        --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name                    accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#            ------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
seq1                               -           500 model_A               PF00001     480     1e-10   50.0   0.1   1   1     2e-10   3e-10   49.0   0.1     5   480     3   478     2   479 0.95 a domtblout description with spaces
"""  # noqa: E501
TSV_DOMTBLOUT_D = """\
#                                                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name                                        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#            ----------------------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
seq4_with_a_much_longer_identifier   -           620 model_D               PF00004     600     2e-15   80.0   0.2   1   1     3e-15   4e-15   79.0   0.2    10   600     8   598     6   601 0.91 a second file description
"""  # noqa: E501
TSV_SPACE_LABEL_A = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 a multi word description
"""  # noqa: E501
TSV_SPACE_LABEL_B = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq2                               -          model_B               PF00002       3e-08   40.2   0.0      4e-08   39.9   0.0   1.0   1   0   0   1   1   1   1 another one
"""  # noqa: E501
TSV_NO_HEADER = """\
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 a multi word description
"""  # noqa: E501
TSV_NHMMER = """\
#target name  accession query name  accession   hmmfrom  hmm to alifrom ali to envfrom  env to  modlen strand   E-value  score  bias description of target
#------------ --------- ------------ --------- ------- ------- ------- ------- ------- ------- ------- ------ --------- ------ ----- ---------------------
seq1          -         model_A      -             1     100       1     100       1     100     100      +     1e-10   50.0   0.1 a description
"""  # noqa: E501
TSV_FOOTER = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 a multi word description
#
# Program:         hmmsearch
# Version:         3.4 (Aug 2023)
# Pipeline mode:   SEARCH
# Query file:      model_A.hmm
# Target file:     all.faa.gz
# Date:            Sun Sep 13 14:01:33 2026
# [ok]
"""  # noqa: E501
TSV_TWO_RUNS = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 first run's hit
#
# Program:         hmmsearch
# Version:         3.4 (Aug 2023)
# Pipeline mode:   SEARCH
# Query file:      model_A.hmm
# Target file:     all.faa.gz
# Date:            Sun Sep 13 14:01:33 2026
# [ok]
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq2                               -          model_B               PF00002       2e-10   49.0   0.1      3e-10   48.0   0.1   1.0   1   0   0   1   1   1   1 second run's hit
#
# Program:         hmmsearch
# [ok]
"""  # noqa: E501
TSV_TRAILING_BLANK = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 a description
#
# Program:         hmmsearch
# [ok]

"""  # noqa: E501
TSV_ZERO_HITS = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
#
# [ok]
"""  # noqa: E501
TSV_HASH_ID = """\
#                                                                           --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                    accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#            ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
seq1                               -          model_A               PF00001       1e-10   50.0   0.1      2e-10   49.0   0.1   1.0   1   0   0   1   1   1   1 first hit
#weird_id                          -          model_A               PF00001       2e-10   49.0   0.1      3e-10   48.0   0.1   1.0   1   0   0   1   1   1   1 hash-prefixed id
#
# Program:         hmmsearch
# [ok]
"""  # noqa: E501


async def table_file(name: str, text: str | bytes) -> File:
    return await File.from_local(local_file(name, text.encode() if isinstance(text, str) else text))


async def tsv_lines(tsv: File) -> list[str]:
    """The TSV's lines as upstream's assertions read them (Java's readLines: a line ends at CR, LF or CRLF)."""
    return gzip.decompress(await file_bytes(tsv)).decode().splitlines()


async def format_synthetic(texts: dict[str, str | bytes], labels: list[str], table_format: str) -> File:
    tables = [await table_file(name, text) for name, text in texts.items()]
    out = await format_tsv(tables, labels, table_format, prefix="test")
    assert out.name == f"test.{table_format}.tsv.gz", out.name
    return out


def rows_labelled(lines: list[str], label: str) -> list[str]:
    return [line for line in lines if line.startswith(label + "\t")]


async def check_real_tables(table_format: str, md5: str, counts: dict[str, int], fields: int) -> None:
    tables = [await fixture(PHYLOPLACE + "hmmrank/" + n) for n in NRDJ_TABLES[table_format]]
    out = await format_tsv(tables, ["NrdJ", "NrdJm1"], table_format, prefix="test")
    assert out.name == f"test.{table_format}.tsv.gz", out.name
    await assert_gunzipped_md5(out, md5, label=out.name)
    lines = await tsv_lines(out)
    for label, count in counts.items():
        assert len(rows_labelled(lines, label)) == count, f"{len(rows_labelled(lines, label))} {label} rows"
    assert all(len(line.split("\t")) == fields for line in lines[1:]), "rows with the wrong field count"


@env.task
async def test_hmmer_format_tsv_real_tblout() -> None:
    await check_real_tables("tblout", "d46c71d6d56f735ee4b23028a79f3013", {"NrdJ": 4, "NrdJm1": 5}, 20)


@env.task
async def test_hmmer_format_tsv_real_domtblout() -> None:
    await check_real_tables("domtblout", "4845163f3747f16d1d0a85782b34e2f9", {"NrdJ": 16, "NrdJm1": 12}, 24)


@env.task
async def test_hmmer_format_tsv_labels_not_filenames() -> None:
    out = await format_synthetic(
        {"profileA.tbl": TSV_TBLOUT_A, "profileB.tbl": TSV_TBLOUT_B}, ["profileA", "profileB"], "tblout"
    )
    await assert_gunzipped_md5(out, "e5f6696f24f306b857ca9f571fdbe6fb", label=out.name)
    lines = await tsv_lines(out)
    assert len(rows_labelled(lines, "profileA")) == 2 and len(rows_labelled(lines, "profileB")) == 1, lines
    assert any("a multi word description" in line for line in lines), lines


@env.task
async def test_hmmer_format_tsv_domtblout_header_widths() -> None:
    out = await format_synthetic(
        {"profileC.domtbl": TSV_DOMTBLOUT_C, "profileD.domtbl": TSV_DOMTBLOUT_D}, ["profileC", "profileD"], "domtblout"
    )
    await assert_gunzipped_md5(out, "3686c09d5c0a932c8add906dc37e9f23", label=out.name)
    lines = await tsv_lines(out)
    assert len(lines[0].split("\t")) == 24, lines[0]
    assert len(rows_labelled(lines, "profileC")) == 1 and len(rows_labelled(lines, "profileD")) == 1, lines
    assert any("a domtblout description with spaces" in line for line in lines), lines
    assert any("a second file description" in line for line in lines), lines
    assert len(rows_labelled(lines, "profileD")[0].split("\t")) == 24, rows_labelled(lines, "profileD")


@env.task
async def test_hmmer_format_tsv_label_with_space() -> None:
    out = await format_synthetic(
        {"profileA.tbl": TSV_SPACE_LABEL_A, "profileB.tbl": TSV_SPACE_LABEL_B}, ["profile A", "profileB"], "tblout"
    )
    lines = await tsv_lines(out)
    assert any(line.startswith("profile A\t") and "a multi word description" in line for line in lines), lines
    assert any(line.startswith("profileB\t") and "another one" in line for line in lines), lines


async def check_format_fails(text: str, message: str) -> None:
    try:
        await format_synthetic({"profileA.tbl": text}, ["profileA"], "tblout")
    except Exception as e:  # the formatter exits 1 with its message on stderr
        print(f"format_tsv failed as expected: {e}")
        if "hmmer/formattsv" in str(e):
            assert message in str(e), str(e)
    else:
        raise AssertionError(f"format_tsv succeeded where it should fail with {message!r}")


@env.task
async def test_hmmer_format_tsv_no_header_fails() -> None:
    await check_format_fails(TSV_NO_HEADER, "column-count header")


@env.task
async def test_hmmer_format_tsv_other_layout_fails() -> None:
    await check_format_fails(TSV_NHMMER, "expected 19 tblout columns")


@env.task
async def test_hmmer_format_tsv_footer_dropped() -> None:
    lines = await tsv_lines(await format_synthetic({"profileA.tbl": TSV_FOOTER}, ["profileA"], "tblout"))
    assert len(lines) == 2, lines
    assert any(line.startswith("profileA\t") and "a multi word description" in line for line in lines), lines


@env.task
async def test_hmmer_format_tsv_two_runs_fails() -> None:
    await check_format_fails(TSV_TWO_RUNS, "more than one HMMER run")


@env.task
async def test_hmmer_format_tsv_trailing_blank_line() -> None:
    out = await format_synthetic({"profileA.tbl": TSV_TRAILING_BLANK}, ["profileA"], "tblout")
    await assert_gunzipped_md5(out, "41c7183ae09bef3f687cd83d2aba7de9", label=out.name)


@env.task
async def test_hmmer_format_tsv_zero_hits() -> None:
    lines = await tsv_lines(await format_synthetic({"profileA.tbl": TSV_ZERO_HITS}, ["profileA"], "tblout"))
    assert len(lines) == 1, lines


@env.task
async def test_hmmer_format_tsv_hash_target_name() -> None:
    lines = await tsv_lines(await format_synthetic({"profileA.tbl": TSV_HASH_ID}, ["profileA"], "tblout"))
    assert len(lines) == 3, lines
    assert any(line.startswith("profileA\t#weird_id\t") for line in lines), lines


# Upstream's CRLF table (the trailing-blank table's header and row): each line + CR, joined with LF, then CRLF.
TSV_CRLF_LINES = TSV_TRAILING_BLANK.splitlines()[:4]


@env.task
async def test_hmmer_format_tsv_crlf() -> None:
    crlf = ("\n".join(line + "\r" for line in TSV_CRLF_LINES) + "\r\n").encode()
    lines = await tsv_lines(await format_synthetic({"profileA.tbl": crlf}, ["profileA"], "tblout"))
    assert len(lines) == 2, lines
    assert all("\r" not in line for line in lines), lines
    assert any(line.startswith("profileA\t") and "a description" in line for line in lines), lines


tests = [
    test_hmmsearch_txt,
    test_hmmsearch_txt_sto_tbl_domtbl,
    test_hmmalign,
    test_hmmbuild,
    test_hmmrank_16s,
    test_hmmrank_16s_domain_coordinates,
    test_hmmrank_known_domain_coordinates,
    test_hmmrank_several_models_in_one_table,
    test_hmmrank_hit_without_domains_is_na,
    test_hmmrank_tied_domain_starts,
    test_hmmrank_prefix_with_quote,
    test_hmmrank_prefix_with_shell_characters,
    test_hmmrank_prefix_with_newline_and_heredoc_delimiter,
    test_hmmer_esl_sfetch_index,
    test_hmmer_esl_alimask_no_masks,
    test_hmmer_esl_alimask_all_masks,
    test_hmmer_esl_reformat_afa,
    test_hmmer_esl_reformat_phylip,
    test_hmmer_esl_reformat_unalign,
    test_hmmer_format_tsv_real_tblout,
    test_hmmer_format_tsv_real_domtblout,
    test_hmmer_format_tsv_labels_not_filenames,
    test_hmmer_format_tsv_domtblout_header_widths,
    test_hmmer_format_tsv_label_with_space,
    test_hmmer_format_tsv_no_header_fails,
    test_hmmer_format_tsv_other_layout_fails,
    test_hmmer_format_tsv_footer_dropped,
    test_hmmer_format_tsv_two_runs_fails,
    test_hmmer_format_tsv_trailing_blank_line,
    test_hmmer_format_tsv_zero_hits,
    test_hmmer_format_tsv_hash_target_name,
    test_hmmer_format_tsv_crlf,
]
