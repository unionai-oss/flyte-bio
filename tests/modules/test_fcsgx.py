"""Tests for flyte_bio.modules.fcsgx.

The cases and expected values are upstream's fcsgx/fetchdb, fcsgx/rungx and
fcsgx/cleangenome tests (modules commit 9a1cc9f, FCS-GX 0.5.5), against NCBI's
4.6 GB test database (``test-only``), as upstream. As upstream, the database
is checked by its file list (too large to md5) and the run_gx.py reports by
the lines that carry no timestamp or path.
"""

import zlib

from flyte.io import File

from flyte_bio.modules.fcsgx import fcsgx_cleangenome, fcsgx_fetchdb, fcsgx_rungx
from tests.framework import assert_md5, env, fixture

TEST_MANIFEST = "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest"
SARSCOV2_TAXID = "2697049"
COMBO_FASTA = "https://zenodo.org/records/10932013/files/FCS_combo_test.fa"
DATABASE_FILES = [
    "test-only.blast_div.tsv.gz",
    "test-only.gxi",
    "test-only.gxs",
    "test-only.manifest",
    "test-only.meta.jsonl",
    "test-only.seq_info.tsv.gz",
    "test-only.taxa.tsv",
]


async def tar_gz_members(archive: File) -> list[str]:
    """The regular files' paths in a ``.tar.gz``, streamed (the database is too large to hold)."""
    inflate = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buffer = bytearray()
    names: list[str] = []
    skip = 0  # payload bytes still to drop

    def consume() -> bool:
        nonlocal skip
        while True:
            if skip:
                dropped = min(skip, len(buffer))
                del buffer[:dropped]
                skip -= dropped
                if skip:
                    return True
            if len(buffer) < 512:
                return True
            header = bytes(buffer[:512])
            del buffer[:512]
            if header == bytes(512):
                return False
            name = header[:100].rstrip(b"\0").decode()
            prefix = header[345:500].rstrip(b"\0").decode()
            size = int(header[124:136].rstrip(b"\0 ").decode() or "0", 8)
            if header[156:157] in (b"0", b"\0"):
                names.append(f"{prefix}/{name}" if prefix else name)
            skip = (size + 511) // 512 * 512

    async with archive.open("rb") as fh:
        while True:
            chunk = await fh.read(1 << 22)
            if not chunk:
                break
            buffer.extend(inflate.decompress(bytes(chunk)))
            if not consume():
                break
    return sorted(names)


async def lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().split("\n")


@env.task
async def test_fcsgx_fetchdb_test_manifest() -> None:
    # upstream case: fcsgx/fetchdb "fcsgx-test-data - test-manifest"
    database = await fcsgx_fetchdb(TEST_MANIFEST)
    names = await tar_gz_members(database)
    expected = [f"gxdb_test-only/{name}" for name in DATABASE_FILES]
    assert names == expected, names


@env.task
async def test_fcsgx_rungx_sarscov2() -> None:
    # upstream case: fcsgx/rungx "sarscov2 - fasta"
    database = await fcsgx_fetchdb(TEST_MANIFEST)
    genome = await fixture("genomics/sarscov2/genome/genome.fasta")
    result = await fcsgx_rungx(genome, SARSCOV2_TAXID, database, prefix="test")
    report = await lines(result.fcsgx_report)
    assert report[1] == "#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name", report[:2]
    taxonomy = await lines(result.taxonomy_report)
    expected_taxonomy = [
        (
            "#seq-id\tseq-len\t(xp,lc,co,n,mt,pt,pm)-len\tcvg-by-all\tsep1\ttax-name-1\ttax-id-1\tdiv-1\tcvg-by-div-1"
            "\tcvg-by-tax-1\tscore-1\tsep2\ttax-id-2\tdiv-2\tcvg-by-div-2\tcvg-by-tax-2\tscore-2\tsep3\ttax-id-3\tdiv-3"
            "\tcvg-by-div-3\tcvg-by-tax-3\tscore-3\tsep4\ttax-id-4\tdiv-4\tcvg-by-div-4\tcvg-by-tax-4\tscore-4\tsep5"
            "\treserved\tresult\tdiv\tdiv_pct_cvg"
        ),
        (
            "MT192765.1\t29829\t0,0,0,0,0,0,0\t0\t|\t\t\t\t\t\t\t|\t\t\t\t\t\t|\t\t\t\t\t\t|\t\t\t\t\t\t|\tn/a"
            "\tlow-coverage\tnone\t0"
        ),
    ]
    assert taxonomy[1:3] == expected_taxonomy, taxonomy[:3]
    log = await lines(result.log)
    expected_log = [
        "",
        "-----------------------------------------------------------------------------",
        "",
        "tax-id    : 2697049",
        "fasta     : genome.fasta",
        "size      : 0.02 MiB",
        "split-fa  : True",
        "BLAST-div : viruses",
        "gx-div    : virs:viruses",
        "w/same-tax: True",
    ]
    assert log[:10] == expected_log, log[:10]
    log_text = "\n".join(log)
    assert "fcs_gx_report.txt contamination summary:" in log_text, "log lacks the contamination summary"
    assert "fcs_gx_report.txt action summary:" in log_text, "log lacks the action summary"
    assert result.hits is None, "run_gx.py reported hits"


@env.task
async def test_fcsgx_cleangenome_combo() -> None:
    # upstream case: fcsgx/cleangenome "sarscov2 - fasta", on run_gx.py's report for the same FASTA
    database = await fcsgx_fetchdb(TEST_MANIFEST)
    combo = await fixture(COMBO_FASTA)
    screened = await fcsgx_rungx(combo, SARSCOV2_TAXID, database, prefix="test")
    result = await fcsgx_cleangenome(combo, screened.fcsgx_report, prefix="test")
    await assert_md5(result.cleaned, "59ba4b4c773f02a0dc90034762aa1ddf", label="test.cleaned.fasta")
    await assert_md5(result.contaminants, "17c11ea7d2d0076be5183b87618ee91e", label="test.contaminants.fasta")


tests = [test_fcsgx_fetchdb_test_manifest, test_fcsgx_rungx_sarscov2, test_fcsgx_cleangenome_combo]
