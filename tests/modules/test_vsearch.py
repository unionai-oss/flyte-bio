"""Tests for flyte_bio.modules.vsearch.

The cases and expected md5s are upstream's vsearch/cluster tests (modules
commit 9a1cc9f, VSEARCH 2.31.0, samtools 1.23.1). nf-test hashes ``.gz``
outputs decompressed, and the BAM as is. The module's 2 CPUs match upstream's
process_low; the BAM's header records VSEARCH's ``--threads``.

The vsearch_sort cases are upstream's vsearch/sort tests (same commit,
VSEARCH 2.31.0 in its VSEARCH-only image), with their config's
``<id>_sorted`` prefix.
"""

from flyte.io import File

from flyte_bio.modules.vsearch import vsearch_cluster, vsearch_sort
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


tests = [
    test_vsearch_cluster_fast,
    test_vsearch_cluster_size_samout,
    test_vsearch_cluster_smallmem,
    test_vsearch_cluster_unoise,
    test_vsearch_cluster_userout,
    test_vsearch_sort_size,
    test_vsearch_sort_length,
]
