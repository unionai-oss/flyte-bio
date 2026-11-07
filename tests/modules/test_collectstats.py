"""Tests for flyte_bio.modules.collectstats.

The cases and expected values are upstream's custom/collectstats tests
(modules commit 9a1cc9f): two samples (``sample1`` paired-end, ``sample2``
single-end), with and without the optional trimming reports, BBDuk logs and
merge table. As upstream, the table is checked by its (decompressed) md5 and
line by line, and the versions file by its md5.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.collectstats import CollectedStats, collect_stats
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

SAMPLES = [("sample1", False), ("sample2", True)]
FEATURES = [
    "CDS",
    "rRNA",
    "Unassigned_NoFeatures",
    "Unassigned_Ambiguity",
    "Unassigned_MultiMapping",
    "Unassigned_Unmapped",
]
# The header's columns after the trimming and BBDuk counts, and the merge table's.
COUNT_COLUMNS = [
    "idxs_n_mapped",
    "idxs_n_unmapped",
    "CDS",
    "Unassigned_Ambiguity",
    "Unassigned_MultiMapping",
    "Unassigned_NoFeatures",
    "Unassigned_Unmapped",
    "rRNA",
]
MERGETAB_COLUMNS = "n_annotated_taxa\tn_annotated_functions"
VERSIONS_MD5 = "b415587aef40e56f5d6886974db8efa3"


async def required_inputs() -> dict[str, list[File]]:
    return {
        "idxstats": [await fixture(f"generic/txt/{s}.idxstats") for s, _ in SAMPLES],
        "featurecounts": [
            await fixture(f"generic/tsv/{s}.{feature}.featureCounts.tsv") for s, _ in SAMPLES for feature in FEATURES
        ],
    }


async def check(result: CollectedStats, md5: str, lines: list[str]) -> None:
    assert result.overall_stats.name == "run.overall_stats.tsv.gz", f"wrote {result.overall_stats.name}"
    await assert_gunzipped_md5(result.overall_stats, md5, label="overall stats")
    async with result.overall_stats.open("rb") as fh:
        got = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    assert got == lines, f"overall stats: {got}"
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


@env.task
async def test_collectstats_all_optional_inputs() -> None:
    # upstream case: custom/collectstats "two samples - all optional inputs provided"
    result = await collect_stats(
        SAMPLES,
        **await required_inputs(),
        trimlogs=[
            await fixture("generic/txt/sample1_1.fastq.gz_trimming_report.txt"),
            await fixture("generic/txt/sample2.fastq.gz_trimming_report.txt"),
        ],
        bblogs=[await fixture(f"generic/txt/{s}.bbduk.log") for s, _ in SAMPLES],
        mergetab=await fixture("generic/tsv/mergetab.tsv"),
        prefix="run",
    )
    await check(
        result,
        "cd5981bb0f52f78c799fd2168e6eef8b",
        [
            "\t".join(["sample\tn_post_trimming\tn_non_contaminated", *COUNT_COLUMNS, MERGETAB_COLUMNS]),
            "sample1\t600\t580\t550\t15\t200\t3\t7\t15\t20\t30\t42\t17",
            "sample2\t150\t150\t180\t5\t40\t1\t2\t8\t5\t10\t30\t12",
        ],
    )


@env.task
async def test_collectstats_optional_inputs_absent() -> None:
    # upstream case: custom/collectstats "two samples - optional inputs absent"
    result = await collect_stats(SAMPLES, **await required_inputs(), prefix="run")
    await check(
        result,
        "d6557299fc03c8fc59045c5ef665379a",
        [
            "\t".join(["sample\tn_post_trimming", *COUNT_COLUMNS]),
            "sample1\tNA\t550\t15\t200\t3\t7\t15\t20\t30",
            "sample2\tNA\t180\t5\t40\t1\t2\t8\t5\t10",
        ],
    )


tests = [
    test_collectstats_all_optional_inputs,
    test_collectstats_optional_inputs_absent,
]
