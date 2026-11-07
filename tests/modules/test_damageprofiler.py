"""Tests for flyte_bio.modules.damageprofiler.

The cases and expected md5s are upstream's damageprofiler tests (modules
commit 9a1cc9f, DamageProfiler 1.1). As upstream, the plots (``.pdf``,
``.svg``) and the log aren't reproducible, so only their presence is checked;
every other file is checked by md5. All three cases give the same md5s.
"""

from flyte.io import Dir, File

from flyte_bio.modules.damageprofiler import damageprofiler
from tests.framework import assert_md5, env, fixture

HUMAN = "genomics/homo_sapiens/"
BAM = HUMAN + "illumina/bam/test.paired_end.markduplicates.sorted.bam"

DATA_MD5S = {
    "3pGtoA_freq.txt": "56c165bbc9197e0b3d5bd3d464f9e235",
    "3p_freq_misincorporations.txt": "62e60259dd7fa2b54cdab6afc1adc9c4",
    "5pCtoT_freq.txt": "38b5d7c87ec2d8422bc5587095711c5f",
    "5p_freq_misincorporations.txt": "05e041ff47e7c6d08ede51f5c81a2464",
    "DNA_comp_genome.txt": "e4d949c6aceddc2c2ef14de9c4aafde5",
    "DNA_composition_sample.txt": "16be9c91a06fa7defee376230b12862f",
    "dmgprof.json": "a011a47edaadb3a07527c329a6492570",
    "editDistance.txt": "7590dac007633ee11025cc99ec9f009b",
    "lgdistribution.txt": "583fd98094de3fb462b133daeaadcf52",
    "misincorporation.txt": "25119c0fcc8839fdbb3ed10922fb07f4",
}
PLOTS = [
    "DamagePlot.pdf",
    "DamagePlot_five_prime.svg",
    "DamagePlot_three_prime.svg",
    "Length_plot.pdf",
    "Length_plot_combined_data.svg",
    "Length_plot_forward_reverse_separated.svg",
    "edit_distance.pdf",
    "edit_distance.svg",
]


async def files_by_path(results: Dir) -> dict[str, File]:
    """Every file under ``results``, keyed by its path relative to it."""
    root = results.path.rstrip("/") + "/"
    return {f.path.removeprefix(root): f async for f in results.walk()}


async def check_profile(files: dict[str, File], expected_names: list[str], subdir: str = "") -> None:
    names = sorted(name.removeprefix(subdir) for name in files if name.startswith(subdir))
    assert names == sorted(expected_names), f"damageprofiler wrote {names}, expected {sorted(expected_names)}"
    for name, md5 in DATA_MD5S.items():
        await assert_md5(files[subdir + name], md5, label=name)


@env.task
async def test_damageprofiler_bam() -> None:
    # upstream case: "homo_sapiens - bam"
    results = await damageprofiler(await fixture(BAM), prefix="test")
    await check_profile(await files_by_path(results), [*DATA_MD5S, *PLOTS, "DamageProfiler.log"])


@env.task
async def test_damageprofiler_bam_fasta_fai() -> None:
    # upstream case: "homo_sapiens - bam - fasta fai"
    results = await damageprofiler(
        await fixture(BAM),
        fasta=await fixture(HUMAN + "genome/genome.fasta"),
        fai=await fixture(HUMAN + "genome/genome.fasta.fai"),
        prefix="test",
    )
    await check_profile(await files_by_path(results), [*DATA_MD5S, *PLOTS, "DamageProfiler.log"])


@env.task
async def test_damageprofiler_bam_species_list() -> None:
    # upstream case: "homo_sapiens - bam - header" (the species list is the genome's header: one reference)
    results = await damageprofiler(
        await fixture(BAM), species_list=await fixture(HUMAN + "genome/genome.header"), prefix="test"
    )
    files = await files_by_path(results)
    subdirs = sorted({name.split("/", 1)[0] for name in files if "/" in name})
    assert len(subdirs) == 1, f"expected one species subdirectory, got {subdirs} in {sorted(files)}"
    await check_profile(files, [*DATA_MD5S, *PLOTS], subdir=subdirs[0] + "/")


tests = [test_damageprofiler_bam, test_damageprofiler_bam_fasta_fai, test_damageprofiler_bam_species_list]
