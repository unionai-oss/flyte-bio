"""Tests for flyte_bio.modules.spades.

The cases and expected md5s are upstream's spades tests (modules commit
9a1cc9f, SPAdes 4.1.0). As upstream's test config, every case runs with
``--rnaviral``; the module's 4 CPUs match upstream's test CPUs. nf-test
computes the md5s over the decompressed assemblies.

The "pe - hmm" case runs on 2 CPUs: there SPAdes 4.1.0 writes upstream's
scaffolds (its gene-cluster sequences) with 1, 2 or 8 threads, but a different
file with 4 (same header, same assembly graph), reproducibly on this cluster.

Upstream's "pe - pacbio" case passes the Nanopore reads in the Nanopore slot,
so it's the same run as "pe - nanopore"; it is ported as written.
"""

import flyte
from flyte.io import File

from flyte_bio.modules.spades import SpadesResult, read_spades_results, spades, spades_cmd
from tests.framework import assert_gunzipped_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
NANOPORE = "genomics/sarscov2/nanopore/fastq/test.fastq.gz"
SISPA = "https://github.com/nf-core/test-datasets/raw/viralrecon/illumina/sispa/"
ARGS = "--rnaviral"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(
    result: SpadesResult,
    scaffolds: str | None,
    contigs: str | None,
    gfa: str,
    warnings: bool = False,
) -> None:
    """Assert the snapshot: these outputs (None: not written), no transcripts or gene clusters."""
    for label, file, md5 in [("scaffolds", result.scaffolds, scaffolds), ("contigs", result.contigs, contigs)]:
        if md5 is None:
            assert file is None, f"spades wrote unexpected {label}"
        else:
            assert file is not None, f"spades wrote no {label}"
            await assert_gunzipped_md5(file, md5, label=label)
    assert result.gfa is not None, "spades wrote no assembly graph"
    await assert_gunzipped_md5(result.gfa, gfa, label="gfa")
    assert result.transcripts is None and result.gene_clusters is None
    assert "SPAdes pipeline finished" in await text(result.log), "spades log lacks 'SPAdes pipeline finished'"
    if warnings:
        assert result.warnings is not None, "spades wrote no warnings log"


async def paired() -> list[File]:
    return [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]


@env.task
async def test_spades_single_end() -> None:
    result = await spades([await fixture(FASTQ + "test_2.fastq.gz")], prefix="test", args=ARGS)
    await check(
        result,
        scaffolds="65ba6a517c152dbe219bf4b5b92bdad7",
        contigs="65ba6a517c152dbe219bf4b5b92bdad7",
        gfa="4a17fce1225be8a9d3aef3333f14e10b",
    )


@env.task
async def test_spades_paired_end() -> None:
    result = await spades(await paired(), prefix="test", args=ARGS)
    await check(
        result,
        scaffolds=None,
        contigs="70e4a5485dd59566b212a199c31c343b",
        gfa="c34fe8688b404619d36251be0ebf563c",
        warnings=True,
    )


NANOPORE_ASSEMBLY = {
    "scaffolds": "7ddaf03740df422a93fcaffbcd7e9679",
    "contigs": "7ddaf03740df422a93fcaffbcd7e9679",
    "gfa": "8efeac4d563caad9cce256a89bc4bb72",
}


@env.task
async def test_spades_paired_end_pacbio() -> None:
    result = await spades(await paired(), nanopore=await fixture(NANOPORE), prefix="test", args=ARGS)
    await check(result, **NANOPORE_ASSEMBLY, warnings=True)


@env.task
async def test_spades_paired_end_nanopore() -> None:
    result = await spades(await paired(), nanopore=await fixture(NANOPORE), prefix="test", args=ARGS)
    await check(result, **NANOPORE_ASSEMBLY, warnings=True)


@env.task
async def test_spades_paired_end_nanopore_yml() -> None:
    yml = await fixture("delete_me/spades/spades_input_yml.yml")
    result = await spades(await paired(), nanopore=await fixture(NANOPORE), dataset_yml=yml, prefix="test", args=ARGS)
    await check(result, **NANOPORE_ASSEMBLY, warnings=True)


@env.task
async def test_spades_paired_end_hmm() -> None:
    reads = [await fixture(SISPA + "SRR11140744_R1.fastq.gz"), await fixture(SISPA + "SRR11140744_R2.fastq.gz")]
    hmm = await fixture("genomics/sarscov2/genome/proteome.hmm.gz")
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    task = spades_cmd.override(
        resources=flyte.Resources(cpu=2, memory="16Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=2")
    )
    results = await task(
        reads_1=reads[0],
        reads_2=reads[1],
        pacbio=None,
        nanopore=None,
        dataset_yml=None,
        hmm=hmm,
        prefix="test",
        args=ARGS,
    )
    result = await read_spades_results(results, "test")
    await check(
        result,
        scaffolds="ce077d5f3380690f8d9a5fe188f82128",
        contigs=None,
        gfa="719116aaf5bc64820d8cade54b5c053b",
    )


tests = [
    test_spades_single_end,
    test_spades_paired_end,
    test_spades_paired_end_pacbio,
    test_spades_paired_end_nanopore,
    test_spades_paired_end_nanopore_yml,
    test_spades_paired_end_hmm,
]
