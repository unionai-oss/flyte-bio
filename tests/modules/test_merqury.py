"""Tests for flyte_bio.modules.merqury.

The cases and expected md5s are upstream's merqury/merqury tests (modules
commit 9a1cc9f, merqury 1.3), with the meryl count and union-sum steps
chained in as upstream sets them up. As upstream, the plots are checked by
name only.
"""

from flyte.io import File

from flyte_bio.modules.merqury import MerquryResult, merqury
from flyte_bio.modules.meryl import meryl_count, meryl_unionsum
from tests.framework import assert_md5, env, fixture

FASTQ = "genomics/homo_sapiens/illumina/fastq/"
GENOME = "genomics/homo_sapiens/genome/genome.fasta"

# Identical in both cases: (output, file name, md5).
SHARED = [
    ("assembly_only_kmers_bed", "genome_only.bed", "b611f22cde0e410a2ca07c1eefd042d3"),
    ("assembly_only_kmers_wig", "genome_only.wig", "19cf44989af72af597ef80d3489b4882"),
    ("dist_hist", "test.dist_only.hist", "e2e6b54b0febef1f0fcf24cd2afd0b7a"),
    ("spectra_cn_hist", "test.genome.spectra-cn.hist", "f4222933758ecce1ec9ab4de7b169c94"),
    ("spectra_asm_hist", "test.spectra-asm.hist", "b5a4ab1485313915d2ce9ce94aa89457"),
    ("assembly_qv", "test.qv", "6e04952bc182221c8b9e242dc3298808"),
    ("scaffold_qv", "test.genome.qv", "c554315aabcc4207c367805cf3090da3"),
    ("read_ploidy", "test.unionsum.hist.ploidy", "3fe5dc377d8c9562980e1b93d01cac94"),
]
ASM_PLOTS = ["test.spectra-asm.fl.png", "test.spectra-asm.ln.png", "test.spectra-asm.st.png"]


def cn_plots(names: list[str]) -> list[str]:
    return [f"{n}.spectra-cn.{kind}.png" for n in names for kind in ("fl", "ln", "st")]


async def read_db() -> File:
    reads = [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]
    return await meryl_unionsum(await meryl_count(reads, kvalue=21, prefix="test"), kvalue=21, prefix="test")


async def check(result: MerquryResult, stats_md5: str, plots: list[str]) -> None:
    for label, name, md5 in [*SHARED, ("stats", "test.completeness.stats", stats_md5)]:
        found = getattr(result, label)
        files = found if isinstance(found, list) else [found]
        assert [f.name for f in files] == [name], f"{label}: {[f.name for f in files]}, expected [{name}]"
        await assert_md5(files[0], md5, label=name)
    names = sorted(f.name for f in result.plots)
    assert names == sorted(plots), f"plots {names}, expected {sorted(plots)}"


@env.task
async def test_merqury_genome() -> None:
    result = await merqury(await read_db(), await fixture(GENOME), prefix="test")
    await check(result, "d923bf7c98342625d022b2c32be8dd3a", [*cn_plots(["test.genome"]), *ASM_PLOTS])


@env.task
async def test_merqury_genome_trio() -> None:
    maternal = await meryl_count([await fixture(FASTQ + "test2_1.fastq.gz")], kvalue=21, prefix="test")
    paternal = await meryl_count([await fixture(FASTQ + "test2_2.fastq.gz")], kvalue=21, prefix="test")
    result = await merqury(
        await read_db(),
        await fixture(GENOME),
        maternal_db=maternal[0],
        paternal_db=paternal[0],
        prefix="test",
    )
    plots = cn_plots(["test.genome", "test.genome.test.test2_1", "test.genome.test.test2_2"])
    await check(result, "a8c7089c212e9dca2e0c47199b3002f9", [*plots, *ASM_PLOTS, "test.hapmers.blob.png"])


tests = [test_merqury_genome, test_merqury_genome_trio]
