"""Tests for flyte_bio.modules.dragonflye.

The cases and expected md5s are upstream's dragonflye tests (modules commit
9a1cc9f, Dragonflye 1.2.1) on the bacterial assembly pipeline's test data,
with their ``ext.args``. As upstream, only the raw contigs and the unpolished
graph are compared by md5; the final contigs (whose order varies) and the log
(which has timestamps) are checked by name. The module's 4 CPUs match
upstream's test CPUs.
"""

from flyte_bio.modules.dragonflye import DragonflyeResult, dragonflye
from tests.framework import assert_md5, env, fixture

DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/bacass/"


async def check(result: DragonflyeResult, assembler: str, raw: str, gfa: str) -> None:
    await assert_md5(result.raw_contigs, raw, label=f"{assembler}.fasta")
    assert result.raw_contigs.path.endswith(f"/{assembler}.fasta"), f"raw contigs {result.raw_contigs.path}"
    assert result.gfa is not None, "dragonflye wrote no unpolished graph"
    assert result.gfa.path.endswith(f"/{assembler}-unpolished.gfa"), f"graph {result.gfa.path}"
    await assert_md5(result.gfa, gfa, label=f"{assembler}-unpolished.gfa")
    names = [f.path.rsplit("/", 1)[-1] for f in result.contigs]
    assert names == ["test.fa", "test.reoriented.fa"], f"contigs {names}"
    assert result.log.path.endswith("/dragonflye.log")


@env.task
async def test_dragonflye_miniasm() -> None:
    # upstream case: "Dragonflye with miniasm"
    reads = await fixture(DATA + "nanopore/subset15000.fq.gz")
    result = await dragonflye(reads, prefix="test", args="--assembler miniasm --gsize 5000000")
    await check(result, "miniasm", "6b8903ba09592df99f43ed05fda488f6", "40ab03a417eafab0cb4ac2c32bd006e1")


@env.task
async def test_dragonflye_raven() -> None:
    # upstream case: "Dragonflye with raven"
    reads = await fixture(DATA + "nanopore/subset15000.fq.gz")
    result = await dragonflye(reads, prefix="test", args="--assembler raven --gsize 5000000")
    await check(result, "raven", "5e4c032b868945c1309dd1d2d793072e", "dea0c9d418c1c53243b9575e50bd14fa")


@env.task
async def test_dragonflye_hybrid() -> None:
    # upstream case: "Dragonflye hybrid" (its args spell the genome size option "-gsize")
    shortreads = [await fixture(DATA + "A1403KPN_R1.fastq.gz"), await fixture(DATA + "A1403KPN_R2.fastq.gz")]
    longreads = await fixture(DATA + "nanopore/A1403KPN.fq.gz")
    result = await dragonflye(longreads, shortreads, prefix="test", args="--assembler miniasm -gsize 5000000")
    await check(result, "miniasm", "b3b05ee807e75ae035ff55d876da607f", "582c59bd7eb9af1fec673d4ecf1772e9")


tests = [test_dragonflye_miniasm, test_dragonflye_raven, test_dragonflye_hybrid]
