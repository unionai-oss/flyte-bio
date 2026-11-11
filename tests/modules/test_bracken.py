"""Tests for flyte_bio.modules.bracken.

The cases and expected md5s are upstream's bracken/bracken tests (modules
commit 9a1cc9f, bracken 3.1, which reports itself as 3.0.1). As upstream,
each case first classifies the reads with kraken2 against the test
Kraken2/Bracken database, then runs Bracken on that report with the test's
``ext.args``: ``-l S -t 10 -r 150``, or ``-l G -t 10 -r 150`` for the genus
case.

The build case is upstream's bracken/build test (same commit): the
SARS-CoV-2 genome and proteome added to a Kraken2 library over the test
protein taxonomy, built with kraken2 without cleaning (and, as that test,
without ``--protein``), then with ``bracken-build``'s defaults. As upstream,
it checks the database directory's listing.

The combine case is upstream's bracken/combinebrackenoutputs test (same
commit): test_1.fastq.gz classified with kraken2 and re-estimated with
Bracken's defaults (upstream's setup channel holds one sample, ``test``),
then combined as ``db``.
"""

import tarfile

from flyte.io import Dir, File

from flyte_bio.modules.bracken import bracken, bracken_build, bracken_combinebrackenoutputs
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.kraken2 import kraken2, kraken2_add, kraken2_build
from tests.framework import assert_md5, env, fixture, fixture_dir

FASTQ = "genomics/sarscov2/illumina/fastq/"
SPECIES_ARGS = "-l S -t 10 -r 150"
GENUS_ARGS = "-l G -t 10 -r 150"


async def database() -> Dir:
    return await fixture_dir("genomics/sarscov2/genome/db/kraken2_bracken.tar.gz", "kraken2_bracken")


async def kraken_report(db: Dir, paired: bool) -> File:
    reads_1 = await fixture(FASTQ + "test_1.fastq.gz")
    reads_2 = await fixture(FASTQ + "test_2.fastq.gz") if paired else None
    return (await kraken2(reads_1, db, reads_2=reads_2, prefix="test")).report


async def check(paired: bool, args: str, tsv_md5: str, report_md5: str) -> None:
    db = await database()
    result = await bracken(await kraken_report(db, paired), db, prefix="test", args=args)
    assert result.report.name == "test.tsv", result.report.name
    await assert_md5(result.report, tsv_md5, label="test.tsv")
    await assert_md5(result.kraken_style_report, report_md5, label="test.kraken2.report_bracken.txt")


@env.task
async def test_bracken_single_end() -> None:
    await check(False, SPECIES_ARGS, "4a21ae14ff8a0311d55f139af5247838", "ca0fbeedc4353b5fdd081688823a33df")


@env.task
async def test_bracken_paired_end_genus() -> None:
    await check(True, GENUS_ARGS, "f609b09d6edb5ebc1ea1435d1dd46cde", "2ce58814420a3690da1f08e10e8d3a30")


@env.task
async def test_bracken_paired_end() -> None:
    await check(True, SPECIES_ARGS, "4a21ae14ff8a0311d55f139af5247838", "ca0fbeedc4353b5fdd081688823a33df")


METAGENOME = "genomics/sarscov2/metagenome/"
BUILT_DATABASE = [
    "database.kraken",
    "database100mers.kmer_distrib",
    "database100mers.kraken",
    "hash.k2d",
    "library",
    "opts.k2d",
    "seqid2taxid.map",
    "taxo.k2d",
    "taxonomy",
]


@env.task
async def test_bracken_build() -> None:
    # upstream case: bracken/build "bracken - db"
    library = await kraken2_add(
        [
            await fixture("genomics/sarscov2/genome/genome.fasta"),
            await fixture("genomics/sarscov2/genome/proteome.fasta"),
        ],
        taxonomy_names=await fixture(METAGENOME + "prot_names.dmp"),
        taxonomy_nodes=await fixture(METAGENOME + "prot_nodes.dmp"),
        accession2taxid=await gunzip(await fixture(METAGENOME + "prot.accession2taxid.gz")),
        prefix="test",
    )
    kraken2_db = await kraken2_build(library, clean=False)
    assert kraken2_db.archive is not None, "kraken2_build(clean=False) wrote no archive"
    result = await bracken_build(kraken2_db.archive)
    with tarfile.open(await result.archive.download()) as tar:
        names = tar.getnames()
    top = sorted({name.split("/")[1] for name in names if name.count("/") >= 1})
    assert top == BUILT_DATABASE, f"bracken-database: {top}"
    distrib = await result.db.get_file("database100mers.kmer_distrib")
    assert distrib is not None, "bracken_build's db has no database100mers.kmer_distrib"


@env.task
async def test_bracken_combinebrackenoutputs() -> None:
    # upstream case: bracken/combinebrackenoutputs "sarscov2 - fastq - fastq"
    db = await database()
    result = await bracken(await kraken_report(db, False), db, prefix="test")
    combined = await bracken_combinebrackenoutputs([result.report], prefix="db")
    assert combined.name == "db.txt", combined.name
    await assert_md5(combined, "fb22ad13f9685d5bf5026542b9b2a6b8", label="db.txt")


tests = [
    test_bracken_single_end,
    test_bracken_paired_end_genus,
    test_bracken_paired_end,
    test_bracken_build,
    test_bracken_combinebrackenoutputs,
]
