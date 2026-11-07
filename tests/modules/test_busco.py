"""Tests for flyte_bio.modules.busco.

The cases and expected md5s are upstream's busco/busco tests (modules commit
9a1cc9f, BUSCO 6.1.0), all run with ``--tar`` as upstream's test config sets.
Upstream lets BUSCO download its lineage in every run; here each lineage is
fetched once by :func:`busco_download` (cached) and passed in as ``lineages``.
"""

import hashlib

from flyte.io import Dir, File

from flyte_bio.modules.busco import BuscoResult, busco, busco_download
from tests.framework import assert_md5, env, fixture

BACTEROIDES_GENOME = "genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz"
PORTIERA_GENOME = "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/genome/genome.fasta"
PORTIERA_PROTEOME = "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/genome/proteome.fasta"
BACTEROIDES_CONTIGS = "genomics/prokaryotes/bacteroides_fragilis/illumina/fasta/test1.contigs.fa.gz"
HUMAN_GENOME = "genomics/homo_sapiens/genome/genome.fasta"
BACTERIA = "bacteria_odb10"
SHORT_SUMMARY_TXT = [
    "BUSCO version",
    "The lineage dataset is",
    "BUSCO was run in mode",
    "Complete BUSCOs",
    "Missing BUSCOs",
    "Dependencies and versions",
]
SHORT_SUMMARY_JSON = ["one_line_summary", "mode", "dataset"]
LOG = ["DEBUG:busco.run_BUSCO", "Results from dataset", "how to cite BUSCO"]
SEQUENCE_TARS = [
    "single_copy_busco_sequences.tar.gz",
    "multi_copy_busco_sequences.tar.gz",
    "fragmented_busco_sequences.tar.gz",
]


def name(path: str) -> str:
    return path.rstrip("/").rsplit("/", 1)[-1]


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def md5(file: File) -> str:
    async with file.open("rb") as fh:
        return hashlib.md5(bytes(await fh.read())).hexdigest()


async def check_contains(file: File | None, wanted: list[str], label: str) -> None:
    assert file is not None, f"no {label}"
    content = await text(file)
    for expected in wanted:
        assert expected in content, f"{label} lacks {expected!r}"


async def check_md5s(files: list[File], md5s: list[str], label: str) -> None:
    assert len(files) == len(md5s), f"{label}: {len(files)} files, expected {len(md5s)}"
    for file, expected in zip(files, md5s, strict=True):
        await assert_md5(file, expected, label=f"{label} {file.path}")


async def check_summaries(result: BuscoResult, count: int) -> None:
    assert len(result.short_summaries_txt) == count, f"{len(result.short_summaries_txt)} txt short summaries"
    assert len(result.short_summaries_json) == count, f"{len(result.short_summaries_json)} json short summaries"
    for summary in result.short_summaries_txt:
        await check_contains(summary, SHORT_SUMMARY_TXT, name(summary.path))
    for summary in result.short_summaries_json:
        await check_contains(summary, SHORT_SUMMARY_JSON, name(summary.path))


async def check_seq_dirs(result: BuscoResult, count: int) -> None:
    assert len(result.seq_dir) == count, f"{len(result.seq_dir)} busco_sequences dirs, expected {count}"
    for seq_dir in result.seq_dir:
        names = [name(f.path) for f in await seq_dir.list_files()]
        for tar in SEQUENCE_TARS:
            assert tar in names, f"{seq_dir.path} lacks {tar} (has {names})"


async def dir_md5s(directory: Dir) -> dict[str, str]:
    return {name(f.path): await md5(f) async for f in directory.walk()}


@env.task
async def test_busco_genome_single_fasta() -> None:
    lineages = await busco_download(BACTERIA)
    result = await busco(
        [await fixture(BACTEROIDES_GENOME)], "genome", BACTERIA, lineages, prefix=f"test-{BACTERIA}", args="--tar"
    )
    await check_summaries(result, 1)
    await assert_md5(result.batch_summary, "12e911830d66bab6dbf3523ac4392597", label="batch_summary")
    await check_md5s(result.full_table, ["327aa7603dca43e7c4adf85f89a85e8e"], "full_table")
    await check_md5s(result.missing_busco_list, ["268d346a92e1297999d96c98a6e7c260"], "missing_busco_list")
    await check_seq_dirs(result, 1)
    await check_contains(result.log, LOG, "log")
    assert result.single_copy_proteins == [] and result.translated_dir == []


@env.task
async def test_busco_genome_multi_fasta() -> None:
    lineages = await busco_download(BACTERIA)
    fasta = [await fixture(BACTEROIDES_GENOME), await fixture(PORTIERA_GENOME)]
    result = await busco(fasta, "genome", BACTERIA, lineages, prefix=f"test-{BACTERIA}", args="--tar")
    await check_summaries(result, 2)
    # Sorted by path: genome.fasta (Portiera) before genome.fna (Bacteroides), as upstream's snapshot.
    await check_md5s(
        result.full_table, ["223cff423706f07a121f53115bcd6d69", "327aa7603dca43e7c4adf85f89a85e8e"], "full_table"
    )
    await check_md5s(
        result.missing_busco_list,
        ["9a8eb41928d1c21677b9d38abaf2e790", "268d346a92e1297999d96c98a6e7c260"],
        "missing_busco_list",
    )
    await check_seq_dirs(result, 2)
    await check_contains(result.log, LOG, "log")
    assert result.single_copy_proteins == [] and result.translated_dir == []


@env.task
async def test_busco_eukaryote_augustus() -> None:
    lineage = "eukaryota_odb10"
    lineages = await busco_download(lineage)
    result = await busco(
        [await fixture(HUMAN_GENOME)], "genome", lineage, lineages, prefix=f"test-{lineage}", args="--tar --augustus"
    )
    await assert_md5(result.batch_summary, "3ea3bdc423a461dae514d816bdc61c89", label="batch_summary")
    await check_contains(result.log, ["DEBUG:busco.run_BUSCO", "Augustus did not recognize any genes"], "log")
    assert result.short_summaries_json == [] and result.short_summaries_txt == []
    assert result.missing_busco_list == [] and result.full_table == []
    assert result.single_copy_proteins == [] and result.translated_dir == []


@env.task
async def test_busco_protein() -> None:
    lineages = await busco_download(BACTERIA)
    result = await busco(
        [await fixture(PORTIERA_PROTEOME)], "proteins", BACTERIA, lineages, prefix=f"test-{BACTERIA}", args="--tar"
    )
    await check_summaries(result, 1)
    await assert_md5(result.batch_summary, "942dbb2d8ff26240860a794213db14a8", label="batch_summary")
    await check_md5s(result.full_table, ["6cea99feda083a56bf8eb6d9bd19f4be"], "full_table")
    await check_md5s(result.missing_busco_list, ["9a8eb41928d1c21677b9d38abaf2e790"], "missing_busco_list")
    await check_seq_dirs(result, 1)
    await check_contains(result.log, LOG, "log")
    assert result.single_copy_proteins == [] and result.translated_dir == []


TRANSLATED_PROTEINS = {
    "1024388at2.faa": "797d603d262a6595a112e25b73e878b0",
    "1054741at2.faa": "cd4b928cba6b19b4437746ba507e7195",
    "1093223at2.faa": "df9549708e5ffcfaee6a74dd70a0e5dc",
    "1151822at2.faa": "12726afc1cdc40c13392e1596e93df3a",
    "143460at2.faa": "d887431fd988a5556a523440f02d9594",
    "1491686at2.faa": "d03362d19979b27306c192f1c74a84e5",
    "1504821at2.faa": "4f5f6e5c57bac0092c1d85ded73d7e67",
    "1574817at2.faa": "1153e55998c2929eacad2aed7d08d248",
    "1592033at2.faa": "bb7a59e5f3a57ba12d10dabf4c77ab57",
    "1623045at2.faa": "8fe38155feb1802beb97ef7714837bf5",
    "1661836at2.faa": "6c6d592c2fbb0d7a4e5e1f47a15644f0",
    "1674344at2.faa": "bb41b44e53565a54cadf0b780532fe08",
    "1698718at2.faa": "f233860000028eb00329aa85236c71e5",
    "1990650at2.faa": "34a2d29c5f8b6253159ddb7a43fa1829",
    "223233at2.faa": "dec6705c7846c989296e73942f953cbc",
    "402899at2.faa": "acc0f271f9a586d2ce1ee41669b22999",
    "505485at2.faa": "aa0391f8fa5d9bd19b30d844d5a99845",
    "665824at2.faa": "47f8ad43b6a6078206feb48c2e552793",
    "776861at2.faa": "f8b90c13f7c6be828dea3bb920195e3d",
    "874197at2.faa": "8d22a35a768debe6f376fc695d233a69",
    "932854at2.faa": "2eff2de1ab83b22f3234a529a44e22bb",
    "95696at2.faa": "247bfd1aef432f7b5456307768e9149c",
}


@env.task
async def test_busco_transcriptome() -> None:
    lineages = await busco_download(BACTERIA)
    result = await busco(
        [await fixture(BACTEROIDES_CONTIGS)],
        "transcriptome",
        BACTERIA,
        lineages,
        prefix=f"test-{BACTERIA}",
        args="--tar",
    )
    await check_summaries(result, 1)
    await assert_md5(result.batch_summary, "8734b3f379c4c0928e5dd4ea1873dc64", label="batch_summary")
    await check_md5s(result.full_table, ["e9ab16f5dfe319d30828186658c5d02a"], "full_table")
    await check_md5s(result.missing_busco_list, ["5d5bbb10940f6cf862c0f4c979980a21"], "missing_busco_list")
    assert len(result.translated_dir) == 1, f"{len(result.translated_dir)} translated_proteins dirs"
    actual = await dir_md5s(result.translated_dir[0])
    assert actual == TRANSLATED_PROTEINS, (
        f"translated_proteins differ: {sorted(set(actual.items()) ^ set(TRANSLATED_PROTEINS.items()))}"
    )
    await check_md5s(result.single_copy_proteins, ["73e2c5d6a9b0f01f2deea3cc5f21b764"], "single_copy_proteins")
    await check_seq_dirs(result, 1)
    await check_contains(result.log, LOG, "log")


@env.task
async def test_busco_cleanup() -> None:
    lineages = await busco_download(BACTERIA)
    result = await busco(
        [await fixture(BACTEROIDES_GENOME)],
        "genome",
        BACTERIA,
        lineages,
        clean_intermediates=True,
        prefix=f"test-{BACTERIA}",
        args="--tar",
    )
    await assert_md5(result.batch_summary, "12e911830d66bab6dbf3523ac4392597", label="batch_summary")
    await check_md5s(result.full_table, ["327aa7603dca43e7c4adf85f89a85e8e"], "full_table")
    await check_md5s(result.missing_busco_list, ["268d346a92e1297999d96c98a6e7c260"], "missing_busco_list")
    await check_contains(result.log, LOG, "log")
    base = result.busco_dir.path.rstrip("/")
    paths = [f.path[len(base) :].lstrip("/").split("/") async for f in result.busco_dir.walk()]
    leftovers = [
        "/".join(p)
        for p in paths
        if (len(p) > 2 and p[1] == "auto_lineage")
        or (len(p) > 3 and p[2] in ("miniprot_output", "hmmer_output", ".bbtools_output"))
    ]
    assert leftovers == [], f"intermediates left after cleanup: {leftovers[:5]}"
    assert result.single_copy_proteins == [] and result.translated_dir == []


tests = [
    test_busco_genome_single_fasta,
    test_busco_genome_multi_fasta,
    test_busco_eukaryote_augustus,
    test_busco_protein,
    test_busco_transcriptome,
    test_busco_cleanup,
]
