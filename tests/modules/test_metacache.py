"""Tests for flyte_bio.modules.metacache.

The cases and expected md5s are upstream's metacache/build tests (modules
commit 9a1cc9f, MetaCache 2.5.0): one FASTA, one gzipped FASTA, two FASTAs,
one FASTQ, and one FASTA over the test taxonomy, without and with a
sequence-to-taxon map.
"""

from flyte.io import File

from flyte_bio.modules.metacache import metacache_build
from tests.framework import assert_md5, env, fixture

SARSCOV2 = "genomics/sarscov2/"
TRANSCRIPTOME = "genomics/homo_sapiens/genome/transcriptome.fasta"


async def check_metacache(
    genomes: list[File],
    meta_md5: str,
    cache_md5: str,
    taxonomy: bool = False,
    seq2taxid: bool = False,
) -> None:
    result = await metacache_build(
        genomes,
        taxonomy=[await fixture(SARSCOV2 + "metagenome/names.dmp"), await fixture(SARSCOV2 + "metagenome/nodes.dmp")]
        if taxonomy
        else None,
        seq2taxid=await fixture(SARSCOV2 + "metagenome/seqid2taxid.map") if seq2taxid else None,
        prefix="test",
    )
    assert [c.name for c in result.caches] == ["test.cache0"], [c.name for c in result.caches]
    await assert_md5(result.meta, meta_md5, label="test.meta")
    await assert_md5(result.caches[0], cache_md5, label="test.cache0")


@env.task
async def test_metacache_build_fasta() -> None:
    # upstream case: metacache/build "one fasta"
    await check_metacache(
        [await fixture(TRANSCRIPTOME)], "44bda0e589b01622a6376e410f829907", "1209d774f809cbde2a630f08d99382ac"
    )


@env.task
async def test_metacache_build_fasta_gz() -> None:
    # upstream case: metacache/build "one fasta.gz"
    await check_metacache(
        [await fixture(SARSCOV2 + "genome/genome.fasta.gz")],
        "4cea57b64d29695e824c2c5ab880d2e7",
        "a4357c9b9512635f95d4507620e36fce",
    )


@env.task
async def test_metacache_build_two_fasta() -> None:
    # upstream case: metacache/build "two fasta"
    await check_metacache(
        [await fixture(SARSCOV2 + "genome/genome.fasta"), await fixture(TRANSCRIPTOME)],
        "073f7fa32f14cac528c8e18d68f65e29",
        "2b0ba3d7c81365d62ea39855dab8645c",
    )


@env.task
async def test_metacache_build_fastq() -> None:
    # upstream case: metacache/build "one fastq"
    await check_metacache(
        [await fixture(SARSCOV2 + "illumina/fastq/test_1.fastq.gz")],
        "d27079a160a7be4b7ced25c1c6f33f60",
        "2c3aba2e87bb82e66d7e83434fff8c50",
    )


@env.task
async def test_metacache_build_taxonomy() -> None:
    # upstream case: metacache/build "one fasta + taxonomy"
    await check_metacache(
        [await fixture(SARSCOV2 + "genome/genome.fasta")],
        "7dddd6667a0c7d7d3531d093fee985b7",
        "a4357c9b9512635f95d4507620e36fce",
        taxonomy=True,
    )


@env.task
async def test_metacache_build_taxonomy_seq2taxid() -> None:
    # upstream case: metacache/build "one fasta + taxonomy + seq2taxid"
    await check_metacache(
        [await fixture(SARSCOV2 + "genome/genome.fasta")],
        "7dddd6667a0c7d7d3531d093fee985b7",
        "a4357c9b9512635f95d4507620e36fce",
        taxonomy=True,
        seq2taxid=True,
    )


tests = [
    test_metacache_build_fasta,
    test_metacache_build_fasta_gz,
    test_metacache_build_two_fasta,
    test_metacache_build_fastq,
    test_metacache_build_taxonomy,
    test_metacache_build_taxonomy_seq2taxid,
]
