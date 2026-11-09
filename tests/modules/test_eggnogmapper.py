"""Tests for flyte_bio.modules.eggnogmapper.

The cases and expected values are upstream's eggnogmapper tests (modules
commit 9a1cc9f, eggNOG-mapper 2.1.13): the SARS-CoV-2 proteome searched with
DIAMOND (against a database :func:`diamond_makedb` builds from it, as
upstream's setup does) and MMseqs2 (a database from :func:`mmseqs_createdb`),
then annotated from the diamond run's seed orthologs (``no_search``). As
upstream, the annotation database is the test-datasets' tiny ``eggnog.db``
(60 KB, annotating one protein), not eggNOG's full one (tens of GiB), alone
in the data directory. Upstream snapshots lines 4-7 of the annotation table
(the first three record the command and date), lines 6-19 of the diamond
run's seed orthologs and its hits' md5, and the mmseqs run's line counts.
"""

import shutil
import tempfile
from pathlib import Path

from flyte.io import Dir

from flyte_bio.modules.diamond import diamond_makedb
from flyte_bio.modules.eggnogmapper import EggnogmapperResult, eggnogmapper
from flyte_bio.modules.mmseqs import mmseqs_createdb
from tests.framework import assert_md5, env, fixture

PROTEOME = "genomics/sarscov2/genome/proteome.fasta"
HEADER = (
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\t"
    "Preferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
    "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs"
)


def annotation_lines(score: str) -> list[str]:
    row = (
        f"ENSSASP00005000002.1\tENSSASP00005000002.1\t0.0\t{score}\t"
        "COG0498@1|root,COG0498@2|Bacteria,1MUWQ@1224|Proteobacteria,2VHR6@28216|Betaproteobacteria,"
        "2KUMA@206389|Rhodocyclales\t1224|Proteobacteria\tE\tthreonine synthase" + "\t-" * 13
    )
    return ["##", HEADER, row, "## 12 queries scanned"]


DIAMOND_ORTHOLOGS = [
    "#qseqid\tsseqid\tevalue\tbitscore\tqstart\tqend\tsstart\tsend\tpident\tqcov\tscov",
    *(
        f"ENSSASP00005{n:06d}.1\tENSSASP00005{n:06d}.1\t{evalue}\t{score}\t1\t{length}\t1\t{length}\t100.0\t100.0\t100.0"
        for n, evalue, score, length in [
            (2, "0.0", "14179.0", 7096),
            (3, "0.0", "8659.0", 4405),
            (4, "0.0", "2469.0", 1273),
            (5, "2.54e-286", "769.0", 419),
            (6, "4.58e-205", "551.0", 275),
            (7, "4.9e-160", "433.0", 222),
            (8, "5.64e-94", "258.0", 121),
            (9, "1.44e-87", "242.0", 121),
            (10, "1.14e-47", "137.0", 75),
            (11, "3.07e-39", "115.0", 61),
            (12, "6.59e-20", "65.5", 43),
            (13, "1.07e-25", "79.7", 38),
        ]
    ),
    "## 12 queries scanned",
]


async def lines(file) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


@env.task(cache="auto")
async def eggnogmapper_test_data_dir() -> Dir:
    """A data directory holding just the test-datasets' ``eggnog.db``, as upstream stages it."""
    db = await (await fixture("delete_me/eggnogmapper/eggnog.db")).download()
    root = Path(tempfile.mkdtemp())
    shutil.copy(db, root / "eggnog.db")
    return await Dir.from_local(str(root))


async def run_diamond() -> EggnogmapperResult:
    proteome = await fixture(PROTEOME)
    db = await diamond_makedb(proteome, prefix="dmnd_db")
    return await eggnogmapper(proteome, "diamond", await eggnogmapper_test_data_dir(), db=db, prefix="test")


async def check_annotations(result: EggnogmapperResult, score: str) -> None:
    assert result.annotations.name == "test.emapper.annotations", f"annotations are {result.annotations.name}"
    found = (await lines(result.annotations))[3:7]
    assert found == annotation_lines(score), f"annotation lines 4-7: {found}"


@env.task
async def test_eggnogmapper_diamond() -> None:
    result = await run_diamond()
    await check_annotations(result, "14179.0")
    assert result.orthologs is not None and result.hits is not None, "no seed orthologs or hits"
    orthologs = (await lines(result.orthologs))[5:19]
    assert orthologs == DIAMOND_ORTHOLOGS, f"seed ortholog lines 6-19: {orthologs}"
    await assert_md5(result.hits, "9ab2e1810b911f2e6c63f8c0f998ee83", label="test.emapper.hits")


@env.task
async def test_eggnogmapper_mmseqs() -> None:
    proteome = await fixture(PROTEOME)
    db = await mmseqs_createdb(proteome, prefix="mmseqs_db")
    result = await eggnogmapper(proteome, "mmseqs", await eggnogmapper_test_data_dir(), db=db, prefix="test")
    await check_annotations(result, "12938.0")
    assert result.orthologs is not None and result.hits is not None, "no seed orthologs or hits"
    orthologs, hits = len(await lines(result.orthologs)), len(await lines(result.hits))
    assert (orthologs, hits) == (24, 12), f"{orthologs} seed ortholog lines and {hits} hit lines, expected 24 and 12"


@env.task
async def test_eggnogmapper_no_search() -> None:
    setup = await run_diamond()
    assert setup.orthologs is not None, "the diamond run wrote no seed orthologs"
    result = await eggnogmapper(
        await fixture(PROTEOME), "no_search", await eggnogmapper_test_data_dir(), db=setup.orthologs, prefix="test"
    )
    await check_annotations(result, "14179.0")


tests = [
    test_eggnogmapper_diamond,
    test_eggnogmapper_mmseqs,
    test_eggnogmapper_no_search,
]
