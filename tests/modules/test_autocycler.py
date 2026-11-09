"""Tests for flyte_bio.modules.autocycler.

The cases and expected md5s are upstream's autocycler/subsample,
autocycler/compress, autocycler/cluster, autocycler/trim, autocycler/resolve
and autocycler/combine tests (modules commit 9a1cc9f, Autocycler 0.5.2).

As upstream, each test runs the chain up to its tool from the SARS-CoV-2
contigs and scaffolds: compress → cluster → trim and resolve per cluster →
combine per cluster. Upstream's snapshots list one output per cluster, keyed
by the cluster's directory (``cluster_001`` … ``cluster_005``); resolve runs
with prefix ``test``. The subsampled reads are compared decompressed, as
nf-test hashes them.
"""

import asyncio

from flyte.io import File

from flyte_bio.modules.autocycler import (
    AutocyclerClusterResult,
    AutocyclerResolveResult,
    AutocyclerTrimResult,
    autocycler_cluster,
    autocycler_combine,
    autocycler_compress,
    autocycler_resolve,
    autocycler_subsample,
    autocycler_trim,
)
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FASTA = "genomics/sarscov2/illumina/fasta/"
CLUSTERS = [f"cluster_00{i}" for i in range(1, 6)]


async def autocycler_test_cluster() -> AutocyclerClusterResult:
    assemblies = [await fixture(FASTA + "contigs.fasta"), await fixture(FASTA + "scaffolds.fasta")]
    compressed = await autocycler_compress(assemblies)
    result = await autocycler_cluster(compressed.gfa)
    assert list(result.clusters) == CLUSTERS, f"clusters {list(result.clusters)}, expected {CLUSTERS}"
    return result


async def autocycler_test_trim() -> list[AutocyclerTrimResult]:
    clustered = await autocycler_test_cluster()
    return list(await asyncio.gather(*(autocycler_trim(gfa) for gfa in clustered.clusters.values())))


async def autocycler_test_resolve() -> list[AutocyclerResolveResult]:
    trimmed = await autocycler_test_trim()
    return list(await asyncio.gather(*(autocycler_resolve(t.gfa, prefix="test") for t in trimmed)))


async def check_each(files: list[File], md5s: list[str], label: str) -> None:
    for name, file, md5 in zip(CLUSTERS, files, md5s, strict=True):
        await assert_md5(file, md5, label=f"{name} {label}")


@env.task
async def test_autocycler_subsample() -> None:
    # upstream case: "Pseudomonas_fluorescens - fastq_gz"; a tiny genome size for enough depth
    reads = await autocycler_subsample(await fixture("genomics/sarscov2/nanopore/fastq/test.fastq.gz"), 300)
    expected = {
        "sample_01.fastq.gz": "7a9ac52d5071893300808fc11091bdbb",
        "sample_02.fastq.gz": "e1e73dc1f31aa8b30834286d13ab500a",
        "sample_03.fastq.gz": "e5cead9b5bbdf541062423065a5b6f51",
        "sample_04.fastq.gz": "312cb2589411da99ecfd50392e999cdc",
    }
    names = [f.path.rsplit("/", 1)[-1] for f in reads]
    assert names == list(expected), f"subsample wrote {names}, expected {list(expected)}"
    for name, file in zip(names, reads):
        await assert_gunzipped_md5(file, expected[name], label=name)


@env.task
async def test_autocycler_compress() -> None:
    # upstream case: "sarscov2 - fasta"
    assemblies = [await fixture(FASTA + "contigs.fasta"), await fixture(FASTA + "scaffolds.fasta")]
    result = await autocycler_compress(assemblies)
    await assert_md5(result.gfa, "e4918bf93dc58755971bec2baf0897e6", label="input_assemblies.gfa")
    await assert_md5(result.stats, "2d4d10d343d16ef57c9dceb12d8c6e1c", label="input_assemblies.yaml")


@env.task
async def test_autocycler_cluster() -> None:
    # upstream case: "sarscov2 - assembly - fasta"
    result = await autocycler_test_cluster()
    await check_each(
        list(result.clusters.values()),
        [
            "73912628f3a2c2d066bc36d14c7229b9",
            "719a197598d0f1499dccbbfe77ece24f",
            "8313fd7f1ab97e957ed1ee4bb20b655d",
            "6296117f3fca264d84ada14f778a23a6",
            "224d8fe4f6350f9e3c5b69c5078e3059",
        ],
        "1_untrimmed.gfa",
    )
    await check_each(
        list(result.cluster_stats.values()),
        [
            "6ac8ce32b804e839499e95279f3945db",
            "44c31adc6f2f2bae83bb3f90d5fe4ab1",
            "1b7a0a5583d9ba65a87546d11019213b",
            "2a1657662c94f7fa22f8c3eac7cd0e42",
            "2ee32f3c9b1bd21d702665f50c86a005",
        ],
        "1_untrimmed.yaml",
    )
    await assert_md5(result.newick, "c0a4ebfeba0c1bb289be9267933a0f37", label="clustering.newick")
    await assert_md5(result.tsv, "30f7a64acefcc55a1863ff306a2901ac", label="clustering.tsv")
    await assert_md5(result.pairwise_distances, "7e266e6d2b0834c802de7909521e197e", label="pairwise_distances.phylip")
    await assert_md5(result.stats, "1e48271631555f4f7a989620d0d4db8c", label="clustering.yaml")


@env.task
async def test_autocycler_trim() -> None:
    # upstream case: "sarscov2 - assemblies - fasta"
    trimmed = await autocycler_test_trim()
    await check_each(
        [t.gfa for t in trimmed],
        [
            "e4339397d2a188c90fc16a10aebba086",
            "b754228734964e5c1a3f006da2de3c37",
            "1bff9016ffb75edc4c73cdeedbc2b7c4",
            "dbbb6784efded7f8cb436fb254c426fd",
            "e08589f6402186b099efeb367add5a09",
        ],
        "2_trimmed.gfa",
    )
    await check_each(
        [t.stats for t in trimmed],
        [
            "236a2cbde4b5c7a3423b65cab2d9204a",
            "a3ce7e8c652fc91da69088cd776c94e6",
            "bf72f5ecfed4757e2114665b3e2ff4e5",
            "b8ca2603d9b0fbec9e0dd53c2d1003dd",
            "de8ad13f8a9b4b911a08a53b9b88e4bd",
        ],
        "2_trimmed.yaml",
    )


@env.task
async def test_autocycler_resolve() -> None:
    # upstream case: "sarscov2 - assemblies - fasta"; bridged, merged and final graphs are identical here
    resolved = await autocycler_test_resolve()
    md5s = [
        "32ab6dfb62c3e0a89e87827d3b216b11",
        "f8ddaec9123e475eab872cfe35767571",
        "f1d19e1dc4d76283358c880d44714bc2",
        "0886b8a7eb8e86875c43f116cd756d7c",
        "f919ba084f3e6e05481843b7acc7470b",
    ]
    await check_each([r.bridged for r in resolved], md5s, "test_3_bridged.gfa")
    await check_each([r.merged for r in resolved], md5s, "test_4_merged.gfa")
    await check_each([r.resolved for r in resolved], md5s, "test_5_final.gfa")


@env.task
async def test_autocycler_combine() -> None:
    # upstream case: "sarscov2 - assemblies -fasta"; as upstream, each cluster is combined on its own
    resolved = await autocycler_test_resolve()
    combined = await asyncio.gather(*(autocycler_combine(r.resolved) for r in resolved))
    await check_each(
        [c.fasta for c in combined],
        [
            "942e9dd63658ae5ee987c3865f64ef97",
            "3b5f0c0e1e480de5682e5956b40013ee",
            "6f06846e831457b1763b0e457b8c101e",
            "a161425fe8e139d75ae54fdd04e56545",
            "8b8d300eb4be010bd97e4a1aebbd61d4",
        ],
        "consensus_assembly.fasta",
    )
    await check_each(
        [c.gfa for c in combined],
        [
            "c69dd262070f9c81a87034f3aaf9888f",
            "6428e5cb960883dcdc4898fa0e972591",
            "30aaa1f58f19ccbf959ac8fa55a78f99",
            "4b946f7cc4d7055fb61bee1925be2bb0",
            "7b2b22be5424893816a4c12f6bc9e1f2",
        ],
        "consensus_assembly.gfa",
    )
    await check_each(
        [c.stats for c in combined],
        [
            "20dfc59b6e81c68d1b74073032cc0d4d",
            "9c2e31d0d42dc7fe715482a19a5768c2",
            "d19724fae1e2b9d0695ec93d25f080bb",
            "d9fcba2ec57942a133c527ae13fc72e6",
            "a935e95b6c41fcee0858902d797e5613",
        ],
        "consensus_assembly.yaml",
    )


tests = [
    test_autocycler_subsample,
    test_autocycler_compress,
    test_autocycler_cluster,
    test_autocycler_trim,
    test_autocycler_resolve,
    test_autocycler_combine,
]
