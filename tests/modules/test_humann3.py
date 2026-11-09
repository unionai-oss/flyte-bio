"""Tests for flyte_bio.modules.humann3.

The cases are upstream's humann3/humann and humann3/regroup tests (modules
commit 9a1cc9f, HUMAnN 3.6.1), against the demo databases on upstream's
funcprofiler test-data branch, passed as their archives. As upstream:

- humann profiles ``test_1.fastq.gz`` with the MetaPhlAn profile of the same
  reads (MetaPhlAn 4.1.1 on its toy database) and checks the pathway abundance
  and coverage tables by md5 (decompressed), and that no reactions table is
  written. The HUMAnN run pins 2 CPUs, as the values were reproduced with.
- regroup regroups ``test_genefamilies.tsv`` into ``uniref90_level4ec`` and
  checks the table by md5.

Upstream tests humann3/renorm only with a stub. Its case here renormalises
``test_genefamilies.tsv`` to copies per million (the default) and checks each
value against the input's, within HUMAnN's 6 significant digits.
"""

import gzip

import flyte
from flyte.io import File

from flyte_bio.modules.humann3 import humann3_humann_cmd, humann3_regroup, humann3_renorm, humann3_result
from flyte_bio.modules.metaphlan import metaphlan_metaphlan
from tests.framework import assert_gunzipped_md5, env, fixture

FUNCPROFILER_DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/funcprofiler/data/database/humann/v3/"
GENEFAMILIES = "genomics/prokaryotes/metagenome/biobakery/test_genefamilies.tsv"
# The cache key ignores resources: the salt keeps these runs apart from default-resource ones.
TWO_CPUS = {"resources": flyte.Resources(cpu=2, memory="8Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=2")}


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_humann3_humann_fastq() -> None:
    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    profile = await metaphlan_metaphlan([reads], await fixture("delete_me/metaphlan4_database.tar.gz"), prefix="test")
    task = humann3_humann_cmd.as_task().override(**TWO_CPUS)
    results = await task(
        input=reads,
        profile=profile.profile,
        nucleotide_db=await fixture(FUNCPROFILER_DATA + "chocophlan_nfDEMO.tar.gz"),
        protein_db=await fixture(FUNCPROFILER_DATA + "uniref_nfDEMO.tar.gz"),
        prefix="test",
        args="",
    )
    result = await humann3_result(results, "test")
    await assert_gunzipped_md5(result.pathabundance, "65f73f474bdb2a7c8b278766ca789a3c", label="pathabundance")
    assert result.pathcoverage is not None, "no pathcoverage"
    await assert_gunzipped_md5(result.pathcoverage, "17c91a2581f34b55179aa4c79c242b1c", label="pathcoverage")
    assert result.reactions is None, "HUMAnN wrote a reactions table"


@env.task
async def test_humann3_regroup_tsv() -> None:
    table = await humann3_regroup(
        await fixture(GENEFAMILIES),
        "uniref90_level4ec",
        await fixture(FUNCPROFILER_DATA + "utility_nfDEMO.tar.gz"),
        prefix="test",
    )
    assert table.name == "test_regroup.tsv.gz", table.name
    await assert_gunzipped_md5(table, "7b00828cbce7a2aa5ece4db335d2fd28", label="regroup")


@env.task
async def test_humann3_renorm_tsv() -> None:
    source = await fixture(GENEFAMILIES)
    table = await humann3_renorm(source, prefix="test")
    assert table.name == "test_renorm.tsv.gz", table.name
    rows = [line.split("\t") for line in (await text(source)).splitlines()]
    async with table.open("rb") as fh:
        renormed = [line.split("\t") for line in gzip.decompress(bytes(await fh.read())).decode().splitlines()]
    assert renormed[0] == rows[0], f"header {renormed[0]}, expected {rows[0]}"
    total = sum(float(row[1]) for row in rows[1:])
    assert [row[0] for row in renormed[1:]] == [row[0] for row in rows[1:]], "features differ"
    for row, out in zip(rows[1:], renormed[1:], strict=True):
        expected = float(row[1]) * 1e6 / total
        assert abs(float(out[1]) - expected) <= expected * 1e-5, f"{row[0]}: {out[1]}, expected about {expected}"


tests = [test_humann3_humann_fastq, test_humann3_regroup_tsv, test_humann3_renorm_tsv]
