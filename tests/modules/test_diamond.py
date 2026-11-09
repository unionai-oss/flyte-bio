"""Tests for flyte_bio.modules.diamond.

The cluster cases are upstream's diamond/cluster tests (modules commit
9a1cc9f, DIAMOND 2.2.2): the human reference proteome, as a FASTA and as a
``.dmnd`` database, with no ``args``. As upstream, they check the cluster
table's line count (one line per protein), not its md5.

The blastp and blastx cases are upstream's diamond/blastp and diamond/blastx
tests at the same commit (DIAMOND 2.2.1), against databases built by
:func:`diamond_makedb` as upstream's setup builds them. They check the hit
tables' md5s (decompressed for ``--compress 1``, as nf-test hashes them),
that DAA output exists, and blastx's log, as upstream does.
"""

from flyte.io import File

from flyte_bio.modules.diamond import diamond_blastp, diamond_blastx, diamond_cluster, diamond_makedb
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, env, fixture

PROTEOME = "proteomics/database/UP000005640_9606"


async def line_count(file: File) -> int:
    async with file.open("rb") as fh:
        return len(bytes(await fh.read()).decode().splitlines())


async def check(db: File) -> None:
    table = await diamond_cluster(db, prefix="human")
    assert table.name == "human.tsv", f"cluster table is {table.name}"
    lines = await line_count(table)
    assert lines == 20610, f"human.tsv has {lines} lines, expected 20610"


@env.task
async def test_diamond_cluster_human_fasta() -> None:
    await check(await fixture(f"{PROTEOME}.fasta"))


@env.task
async def test_diamond_cluster_human_dmnd() -> None:
    await check(await fixture(f"{PROTEOME}.dmnd"))


SARSCOV2 = "genomics/sarscov2/genome/"
TXT_MD5 = "8131b1afd717f3d5f2f2417c5b562e6e"


async def proteome_db() -> File:
    # Upstream's FIRST_DB: makedb of the gzipped proteome, id db1.
    return await diamond_makedb(await fixture(SARSCOV2 + "proteome.fasta.gz"), prefix="db1")


def assert_name(file: File, expected: str) -> None:
    assert file.name == expected, f"output is {file.name}, expected {expected}"


@env.task
async def test_diamond_blastp_txt() -> None:
    db = await proteome_db()
    out = await diamond_blastp(
        await fixture(SARSCOV2 + "proteome.fasta.gz"), db, prefix="txt.db1", outfmt=6, blast_columns="qseqid qlen"
    )
    assert_name(out, "txt.db1.txt")
    await assert_md5(out, TXT_MD5, label="txt.db1.txt")


@env.task
async def test_diamond_blastp_gz_txt() -> None:
    db = await proteome_db()
    out = await diamond_blastp(
        await fixture(SARSCOV2 + "proteome.fasta.gz"), db, prefix="txt.gz.db1", outfmt=6, blast_columns="qseqid qlen"
    )
    assert_name(out, "txt.gz.db1.txt")
    await assert_md5(out, TXT_MD5, label="txt.gz.db1.txt")


@env.task
async def test_diamond_blastp_daa() -> None:
    db = await proteome_db()
    out = await diamond_blastp(await fixture(SARSCOV2 + "proteome.fasta"), db, prefix="daa.db1", outfmt=100)
    assert_name(out, "daa.db1.daa")
    await assert_nonempty(out, label="daa.db1.daa")


@env.task
async def test_diamond_blastp_txt_compressed() -> None:
    # upstream's "txt - gz" case, with its test config's --compress 1
    db = await proteome_db()
    out = await diamond_blastp(
        await fixture(SARSCOV2 + "proteome.fasta"),
        db,
        prefix="txt.gz.db1",
        outfmt=6,
        blast_columns="qseqid qlen",
        args="--compress 1",
    )
    assert_name(out, "txt.gz.db1.txt.gz")
    await assert_gunzipped_md5(out, TXT_MD5, label="txt.gz.db1.txt.gz")


@env.task
async def test_diamond_blastp_two_fastas_two_dbs() -> None:
    db1 = await proteome_db()
    db2 = await diamond_makedb(await fixture(SARSCOV2 + "proteome_test.faa"), prefix="db2")
    fasta1 = await fixture(SARSCOV2 + "proteome.fasta.gz")
    fasta2 = await fixture(SARSCOV2 + "proteome_test.faa")
    expected = {
        ("fasta1", "db1"): "2515cf88590afa32356497e79a51fce9",
        ("fasta1", "db2"): "8b41a8752379f7bd8722258962c598a4",
        ("fasta2", "db1"): "3be947b0e6c69c59491f817b03b1256f",
        ("fasta2", "db2"): "f2461bcecbf4f87cefe952179505aa8f",
    }
    fastas = {"fasta1": fasta1, "fasta2": fasta2}
    dbs = {"db1": db1, "db2": db2}
    for (fasta, db), md5 in expected.items():
        out = await diamond_blastp(fastas[fasta], dbs[db], prefix=f"{fasta}.{db}", outfmt=6)
        assert_name(out, f"{fasta}.{db}.txt")
        await assert_md5(out, md5, label=f"{fasta}.{db}.txt")


async def check_blastx_log(log: File) -> None:
    assert_name(log, "test.log")
    async with log.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert "11 queries aligned." in lines, "blastx log lacks '11 queries aligned.'"


async def transcriptome_db() -> File:
    return await diamond_makedb(await fixture(SARSCOV2 + "proteome.fasta"), prefix="test2")


@env.task
async def test_diamond_blastx_default_txt() -> None:
    # A nonsense extension falls back to tabular output, as upstream checks.
    result = await diamond_blastx(
        await fixture(SARSCOV2 + "transcriptome.fasta"),
        await transcriptome_db(),
        prefix="test",
        out_ext="tfdfdt",
        blast_columns="qseqid qlen",
    )
    assert_name(result.hits, "test.txt")
    await assert_md5(result.hits, "33dc682dabfa44c7089abbc8fe8b84e4", label="test.txt")
    await check_blastx_log(result.log)


@env.task
async def test_diamond_blastx_daa() -> None:
    result = await diamond_blastx(
        await fixture(SARSCOV2 + "transcriptome.fasta"), await transcriptome_db(), prefix="test", out_ext="daa"
    )
    assert_name(result.hits, "test.daa")
    await assert_nonempty(result.hits, label="test.daa")
    await check_blastx_log(result.log)


tests = [
    test_diamond_cluster_human_fasta,
    test_diamond_cluster_human_dmnd,
    test_diamond_blastp_txt,
    test_diamond_blastp_gz_txt,
    test_diamond_blastp_daa,
    test_diamond_blastp_txt_compressed,
    test_diamond_blastp_two_fastas_two_dbs,
    test_diamond_blastx_default_txt,
    test_diamond_blastx_daa,
]
