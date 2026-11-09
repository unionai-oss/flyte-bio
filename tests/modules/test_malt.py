"""Tests for flyte_bio.modules.malt.

The cases and expected md5s are upstream's malt/build tests (modules commit
9a1cc9f, MALT 0.6.2), with the test's ``-J-Xmx6G --sequenceType DNA``: two
SARS-CoV-2 genomes indexed with a MEGAN mapping database (``-mdb``) and with
an accession-to-taxon binary (``-a2t``). The mapping files are large zips
(0.8 and 2.2 GB): they're downloaded with aria2 and extracted with unzip, as
upstream's setup does, on enough disk. As upstream, the log's last line must
report the peak memory, and the index files are checked by md5.
"""

import flyte
from flyte.io import File

from flyte_bio.modules.aria2 import aria2_cmd
from flyte_bio.modules.malt import MaltBuildResult, malt_build
from flyte_bio.modules.unzip import unzip_cmd
from tests.framework import assert_md5, env, fixture

ARGS = "-J-Xmx6G --sequenceType DNA"
MDB_URL = "https://nf-core-test-datasets.s3.amazonaws.com/modules/data/delete_me/malt/megan-nucl-Feb2022.db.zip"
ABIN_URL = (
    "https://ngi-igenomes.s3.eu-west-1.amazonaws.com/test-data/createtaxdb/taxonomy/nucl_acc2tax-Jul2019.abin.zip"
)
SHARED_MD5S = {
    "index0.idx": "0ba3d8bfb7ef28d08e2a005dd3405c55",
    "ref.inf": "042712533a0187b6566db67c6503a71e",
    "taxonomy.idx": "e7ce35e6238f39fa0c236fcf991546e4",
    "taxonomy.map": "5bb3f2192e925bca2e61e4b54f1671e0",
    "taxonomy.tre": "f76fb2d5aa9b0d637234d48175841e0e",
}
BIG_DISK = flyte.Resources(cpu=1, memory="4Gi", disk="40Gi")


async def malt_mapping_file(url: str, suffix: str) -> File:
    """Download and extract a zipped mapping file; returns the file ending in ``suffix``."""
    name = url.rsplit("/", maxsplit=1)[-1]
    download = aria2_cmd.as_task().override(resources=BIG_DISK)
    downloaded = await (await download(source_url=url, file_name=name, args="-x 4")).get_file(name)
    assert downloaded is not None, f"aria2 downloaded no {name}"
    extract = unzip_cmd.as_task().override(resources=BIG_DISK)
    extracted = await extract(archive=downloaded, args="")
    files = [f async for f in extracted.walk(recursive=True) if f.name.endswith(suffix)]
    assert len(files) == 1, f"{name}: expected one *{suffix}, got {[f.name for f in files]}"
    return files[0]


async def malt_test_build(mapping: File, map_type: str) -> MaltBuildResult:
    genomes = [
        await fixture("genomics/prokaryotes/metagenome/fasta/sarscov2.fasta"),
        await fixture("genomics/sarscov2/genome/genome.fasta"),
    ]
    return await malt_build(genomes, mapping_db=mapping, map_type=map_type, args=ARGS)


async def check_malt(result: MaltBuildResult, expected: dict[str, str]) -> None:
    async with result.log.open("rb") as fh:
        last = bytes(await fh.read()).decode().rstrip("\n").splitlines()[-1]
    assert "Peak memory" in last, f"malt-build.log's last line: {last!r}"
    for name, md5 in expected.items():
        file = await result.index.get_file(name)
        assert file is not None, f"malt_index has no {name}"
        await assert_md5(file, md5, label=name)


@env.task
async def test_malt_build_mdb() -> None:
    # upstream case: malt/build "sarscov2 - fastq - mdb"
    result = await malt_test_build(await malt_mapping_file(MDB_URL, ".db"), "mdb")
    await check_malt(
        result,
        SHARED_MD5S | {"ref.db": "6b36ae031c49feaae50f4cea07d9c7f4", "ref.idx": "8ba66cdf65181c7efee1d366574cb9d7"},
    )


@env.task
async def test_malt_build_abin() -> None:
    # upstream case: malt/build "sarscov2 - fastq - abin"
    result = await malt_test_build(await malt_mapping_file(ABIN_URL, ".abin"), "a2t")
    await check_malt(
        result,
        SHARED_MD5S | {"ref.db": "6b0b40df6094a735fa2843e0ef916a02", "ref.idx": "42b0bc56c1562594c7500f2f86a3ca71"},
    )


tests = [test_malt_build_mdb, test_malt_build_abin]
