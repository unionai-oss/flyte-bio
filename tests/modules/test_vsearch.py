"""Tests for flyte_bio.modules.vsearch.

The cases and expected md5s are upstream's vsearch/cluster tests (modules
commit 9a1cc9f, VSEARCH 2.31.0, samtools 1.23.1). nf-test hashes ``.gz``
outputs decompressed, and the BAM as is. The module's 2 CPUs match upstream's
process_low; the BAM's header records VSEARCH's ``--threads``.

The vsearch_sort cases are upstream's vsearch/sort tests (same commit,
VSEARCH 2.31.0 in its VSEARCH-only image), with their config's
``<id>_sorted`` prefix.

The vsearch_sintax and vsearch_usearch_global cases are upstream's
vsearch/sintax and vsearch/usearchglobal tests (same commit and image): the
SARS-CoV-2 genome classified against its transcriptome by md5, and the
transcriptome aligned to the genome at 98.5% identity. Upstream's first
usearch_global case asks for an unknown output ("xcfert"), so falls back to
the alignment, which upstream checks by name and one query's line; the userout
case's lines are compared sorted, as upstream (VSEARCH's output order depends
on its threads).
"""

from flyte.io import File

from flyte_bio.modules.vsearch import vsearch_cluster, vsearch_sintax, vsearch_sort, vsearch_usearch_global
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

READS = "genomics/sarscov2/illumina/fastq/test_1.fastq.gz"


async def check_single(result, name: str) -> File:
    assert result.output is not None and result.clusters == [], "expected one output"
    assert result.output.name == name, f"output named {result.output.name}, expected {name}"
    return result.output


@env.task
async def test_vsearch_cluster_fast() -> None:
    result = await vsearch_cluster(
        await fixture(READS), "--cluster_fast", "--centroids", prefix="test", args="--id 0.8"
    )
    out = await check_single(result, "test.centroids.fasta.gz")
    await assert_gunzipped_md5(out, "7ce841a7bba1c14ced4e1322d73a2d1f", label="cluster_fast centroids")


@env.task
async def test_vsearch_cluster_size_samout() -> None:
    result = await vsearch_cluster(await fixture(READS), "--cluster_size", "--samout", prefix="test", args="--id 0.8")
    out = await check_single(result, "test.bam")
    await assert_md5(out, "51ce92fe4cb652b847e1e9028d955add", label="cluster_size bam")


@env.task
async def test_vsearch_cluster_smallmem() -> None:
    result = await vsearch_cluster(
        await fixture(READS), "--cluster_smallmem", "--centroids", prefix="test", args="--id 0.8 --usersort"
    )
    out = await check_single(result, "test.centroids.fasta.gz")
    await assert_gunzipped_md5(out, "f8a6497f29c55465df843da5c5af39ca", label="cluster_smallmem centroids")


@env.task
async def test_vsearch_cluster_unoise() -> None:
    result = await vsearch_cluster(
        await fixture(READS), "--cluster_unoise", "--centroids", prefix="test", args="--id 0.8 --minsize 1"
    )
    out = await check_single(result, "test.centroids.fasta.gz")
    await assert_gunzipped_md5(out, "e84ed4c7d22e1b8a05400c3751c29efc", label="cluster_unoise centroids")


@env.task
async def test_vsearch_cluster_userout() -> None:
    result = await vsearch_cluster(
        await fixture(READS), "--cluster_fast", "--userout", prefix="test", args="--id 0.8 --userfields query+target+id"
    )
    out = await check_single(result, "test.out.tsv.gz")
    await assert_gunzipped_md5(out, "fb3256dd2327bf3dd4b9695394f5d5a8", label="userout")


async def check_sort(sort: str, md5: str) -> None:
    out = await vsearch_sort(await fixture(READS), sort, prefix="test_sorted")
    assert out.name == "test_sorted.fasta", f"output named {out.name}, expected test_sorted.fasta"
    await assert_md5(out, md5, label=f"{sort} fasta")


@env.task
async def test_vsearch_sort_size() -> None:
    await check_sort("--sortbysize", "86457b12184a7215734ff1e4d868ac51")


@env.task
async def test_vsearch_sort_length() -> None:
    await check_sort("--sortbylength", "e733a34d0c655966408e773de876e5fb")


GENOME = "genomics/sarscov2/genome/genome.fasta"
TRANSCRIPTOME = "genomics/sarscov2/genome/transcriptome.fasta"


async def text_of(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_vsearch_sintax() -> None:
    out = await vsearch_sintax(await fixture(GENOME), await fixture(TRANSCRIPTOME), prefix="test")
    assert out.name == "test.tsv", f"output named {out.name}, expected test.tsv"
    await assert_md5(out, "93bc75fb89343e23fbae971df157b14a", label="test.tsv")


@env.task
async def test_vsearch_usearch_global_unknown_output() -> None:
    out = await vsearch_usearch_global(
        await fixture(TRANSCRIPTOME), await fixture(GENOME), 0.985, "xcfert", prefix="test"
    )
    assert out.name == "test.aln", f"output named {out.name}, expected test.aln"
    assert "Query >lcl|MT192765.1_cds_QIK50427.1_2" in await text_of(out), "alignment lacks the second CDS"


# Upstream's snapshot: the rows, sorted.
USEARCH_GLOBAL_ROWS = [f"lcl|MT192765.1_cds_QIK50{426 + i}.1_{i + 1}\tMT192765.1\t100.0" for i in range(11)]


@env.task
async def test_vsearch_usearch_global_userout() -> None:
    out = await vsearch_usearch_global(
        await fixture(TRANSCRIPTOME), await fixture(GENOME), 0.985, "userout", "query+target+id", prefix="test"
    )
    assert out.name == "test.tsv", f"output named {out.name}, expected test.tsv"
    rows = sorted((await text_of(out)).splitlines())
    assert rows == USEARCH_GLOBAL_ROWS, rows


tests = [
    test_vsearch_cluster_fast,
    test_vsearch_cluster_size_samout,
    test_vsearch_cluster_smallmem,
    test_vsearch_cluster_unoise,
    test_vsearch_cluster_userout,
    test_vsearch_sort_size,
    test_vsearch_sort_length,
    test_vsearch_sintax,
    test_vsearch_usearch_global_unknown_output,
    test_vsearch_usearch_global_userout,
]
