"""Tests for flyte_bio.modules.catpack.

The cases and expected md5s are upstream's catpack tests (modules commit
9a1cc9f, CAT_pack 6.0). As upstream, every case first builds a database with
``catpack_prepare`` (``--common_prefix database``) from the SARS-CoV-2
proteome and its taxonomy, then chains the tools: contigs → add_names →
summarise, and bins → add_names → summarise. As upstream's tests, the
taxonomy files go in swapped (``prot_nodes.dmp`` as ``names``,
``prot_names.dmp`` as ``nodes``); the snapshot's ``names.dmp`` and
``nodes.dmp`` record that. The classification tools' 4 CPUs match upstream's
test CPUs.

Upstream asserts whether the log's last line says "CAT is done!": it does
for contigs, not for bins. ``catpack_download`` has only a stub test upstream
(the databases are many GB), so it isn't tested here.
"""

from flyte.io import File

from flyte_bio.modules.catpack import (
    CatpackClassifyResult,
    CatpackDatabase,
    catpack_addnames,
    catpack_bins,
    catpack_contigs,
    catpack_prepare,
    catpack_summarise,
)
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/"
METAGENOME = "genomics/sarscov2/metagenome/"

ORF2LCA_CONTIGS = "a623b47f20751db12ce18c9ca7ac2536"
CONTIG2CLASSIFICATION = "3c3c79045bf6ae8b1292ae9afa2ec4af"
DIAMOND = "9e2f9c188b183c18dd9572395a48a066"
FAA = "1f8550f87d044d117422ca02827e4d18"
ORF2LCA_BINS = "27e3ca35dc7b977653b5bbf18076fc26"
BIN2CLASSIFICATION = "8cb8c364c1dea229f68e8c7a1747205d"
BINS_GFF = "cb63331a0282175669107585cf4a66c1"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def sarscov2_database() -> CatpackDatabase:
    """Upstream's setup: the SARS-CoV-2 database, with the taxonomy files swapped as upstream's tests pass them."""
    return await catpack_prepare(
        db_fasta=await fixture(GENOME + "proteome.fasta"),
        names=await fixture(METAGENOME + "prot_nodes.dmp"),
        nodes=await fixture(METAGENOME + "prot_names.dmp"),
        acc2tax=await fixture(METAGENOME + "prot.accession2taxid.gz"),
        prefix="test",
        args="--common_prefix database",
    )


async def check_done(result: CatpackClassifyResult, done: bool) -> None:
    lines = (await text(result.log)).splitlines()
    assert ("CAT is done!" in lines[-1]) == done, f"log's last line {lines[-1]!r}: expected done={done}"


async def check_optional(file: File | None, md5: str, label: str) -> None:
    assert file is not None, f"CAT_pack wrote no {label}"
    await assert_md5(file, md5, label=label)


async def check_contigs(result: CatpackClassifyResult) -> None:
    await assert_md5(result.orf2lca, ORF2LCA_CONTIGS, label="ORF2LCA")
    await assert_md5(result.classification, CONTIG2CLASSIFICATION, label="contig2classification")
    await check_optional(result.diamond, DIAMOND, "diamond")
    await check_optional(result.faa, FAA, "faa")
    await check_optional(result.gff, "4fc7a311726723ce5b17cede1dd1059c", "gff")
    await check_done(result, True)


async def check_bins(result: CatpackClassifyResult, premade: bool = False) -> None:
    await assert_md5(result.orf2lca, ORF2LCA_BINS, label="ORF2LCA")
    await assert_md5(result.classification, BIN2CLASSIFICATION, label="bin2classification")
    await check_optional(result.diamond, DIAMOND, "diamond")
    if premade:
        assert result.gff is None, "CAT_pack wrote genes from premade proteins"
    else:
        await check_optional(result.faa, FAA, "faa")
        await check_optional(result.gff, BINS_GFF, "gff")
    await check_done(result, False)


async def genome_contigs(db: CatpackDatabase) -> CatpackClassifyResult:
    return await catpack_contigs(await fixture(GENOME + "genome.fasta"), db, prefix="test")


async def genome_bins(db: CatpackDatabase) -> CatpackClassifyResult:
    return await catpack_bins(
        await fixture(GENOME + "genome.fasta"), db, bin_suffix=".fasta", prefix="test", args="--bin_suffix .fasta"
    )


@env.task
async def test_catpack_prepare() -> None:
    db = await sarscov2_database()
    names = sorted(f.name for f in await db.db.list_files())
    expected = ["database.dmnd", "database.fastaid2LCAtaxid", "database.taxids_with_multiple_offspring"]
    assert names == expected, f"database files {names}, expected {expected}"
    await assert_md5(await db.taxonomy.get_file("names.dmp"), "c471c27a4ce85ae74d2c63633c9ce1e3", label="names.dmp")
    await assert_md5(await db.taxonomy.get_file("nodes.dmp"), "130f9132095562e09c732679c562f5e9", label="nodes.dmp")


@env.task
async def test_catpack_contigs_fasta() -> None:
    await check_contigs(await genome_contigs(await sarscov2_database()))


@env.task
async def test_catpack_contigs_fasta_gz() -> None:
    db = await sarscov2_database()
    await check_contigs(await catpack_contigs(await fixture(GENOME + "genome.fasta.gz"), db, prefix="test"))


@env.task
async def test_catpack_bins_fasta() -> None:
    await check_bins(await genome_bins(await sarscov2_database()))


@env.task
async def test_catpack_bins_premade_proteins() -> None:
    db = await sarscov2_database()
    contigs = await genome_contigs(db)
    result = await catpack_bins(
        await fixture(GENOME + "genome.fasta"),
        db,
        bin_suffix=".fasta",
        proteins=contigs.faa,
        prefix="test",
        args="--bin_suffix .fasta",
    )
    await check_bins(result, premade=True)


@env.task
async def test_catpack_bins_fasta_gz() -> None:
    db = await sarscov2_database()
    result = await catpack_bins(
        await fixture(GENOME + "genome.fasta.gz"),
        db,
        bin_suffix=".fasta.gz",
        prefix="test",
        args="--bin_suffix .fasta",
    )
    await check_bins(result)


@env.task
async def test_catpack_addnames() -> None:
    db = await sarscov2_database()
    contigs = await genome_contigs(db)
    named = await catpack_addnames(contigs.classification, db.taxonomy, prefix="test")
    await assert_md5(named, "bf9914b7772f148d4ed2604a5eba680f", label="add_names")


@env.task
async def test_catpack_summarise_contigs_with_names() -> None:
    db = await sarscov2_database()
    contigs = await genome_contigs(db)
    named = await catpack_addnames(contigs.classification, db.taxonomy, prefix="test", args="--only_official")
    summary = await catpack_summarise(named, await fixture(GENOME + "genome.fasta"), prefix="test_summary")
    await assert_md5(summary, "f4ae1e889f762d6123e4b8d9e57aa4cc", label="summary")


@env.task
async def test_catpack_summarise_bins_raw() -> None:
    db = await sarscov2_database()
    bins = await genome_bins(db)
    named = await catpack_addnames(bins.classification, db.taxonomy, prefix="test", args="--only_official")
    summary = await catpack_summarise(named, prefix="test_summary")
    await assert_md5(summary, "2fc7ecc57cfe0c3362fba27e9693d266", label="summary")


tests = [
    test_catpack_prepare,
    test_catpack_contigs_fasta,
    test_catpack_contigs_fasta_gz,
    test_catpack_bins_fasta,
    test_catpack_bins_premade_proteins,
    test_catpack_bins_fasta_gz,
    test_catpack_addnames,
    test_catpack_summarise_contigs_with_names,
    test_catpack_summarise_bins_raw,
]
