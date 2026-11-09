"""Tests for flyte_bio.modules.metaphlan.

The cases are upstream's metaphlan/metaphlan tests (modules commit 9a1cc9f,
MetaPhlAn 4.1.1), against upstream's toy database
(``delete_me/metaphlan4_database.tar.gz``), passed as the archive. As
upstream, they check lines 3-6 of the profile, that the BIOM table names its
format, and the bowtie2 mapping by md5; the SAM case checks that a SAM was
written. The runs pin 2 CPUs, as the snapshot values were reproduced with.
"""

import flyte
from flyte.io import File

from flyte_bio.modules.metaphlan import MetaphlanResult, metaphlan_metaphlan_cmd, metaphlan_result
from tests.framework import assert_md5, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
DATABASE = "delete_me/metaphlan4_database.tar.gz"
PROFILE_HEADER = [
    "#SampleID\tMetaphlan_Analysis",
    "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species",
    "UNCLASSIFIED\t-1\t100.0\t",
]
# The cache key ignores resources: the salt keeps these runs apart from default-resource ones.
TWO_CPUS = {"resources": flyte.Resources(cpu=2, memory="8Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=2")}


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def metaphlan_two_cpus(input: list[File], save_samfile: bool = False) -> MetaphlanResult:
    task = metaphlan_metaphlan_cmd.as_task().override(**TWO_CPUS)
    database = await fixture(DATABASE)
    results = await task(input=input, database=database, save_samfile=save_samfile, prefix="test", args="")
    return await metaphlan_result(results, "test")


async def check(result: MetaphlanResult, reads_line: str, bowtie2out_md5: str) -> None:
    lines = (await text(result.profile)).splitlines()[2:6]
    assert lines == [reads_line, *PROFILE_HEADER], f"profile lines 3-6: {lines}"
    biom = (await text(result.biom)).splitlines()
    assert "Biological Observation Matrix" in biom[-1], "the BIOM table doesn't name its format"
    assert result.bowtie2out is not None, "no bowtie2out"
    await assert_md5(result.bowtie2out, bowtie2out_md5, label="bowtie2out")


async def paired() -> list[File]:
    return [await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")]


@env.task
async def test_metaphlan_single_end_fastq() -> None:
    result = await metaphlan_two_cpus([await fixture(FASTQ + "test_1.fastq.gz")])
    await check(result, "#100 reads processed", "ef46a9c6a8ce9cae26fbfd5527116fd5")
    assert result.sam is None


@env.task
async def test_metaphlan_paired_end_fastq() -> None:
    result = await metaphlan_two_cpus(await paired())
    await check(result, "#196 reads processed", "8c1bc21e1d8484b5551bf46331d61bd8")


@env.task
async def test_metaphlan_single_end_fasta() -> None:
    result = await metaphlan_two_cpus([await fixture("genomics/sarscov2/illumina/fasta/contigs.fasta")])
    await check(result, "#4 reads processed", "d99d05d6b011c647adf816f10fc6acbe")


@env.task
async def test_metaphlan_paired_end_fastq_save_sam() -> None:
    result = await metaphlan_two_cpus(await paired(), save_samfile=True)
    await check(result, "#196 reads processed", "8c1bc21e1d8484b5551bf46331d61bd8")
    assert result.sam is not None, "no SAM"
    first = (await text(result.sam)).splitlines()[0]
    assert first.startswith("@"), f"the SAM doesn't start with a header: {first!r}"


tests = [
    test_metaphlan_single_end_fastq,
    test_metaphlan_paired_end_fastq,
    test_metaphlan_single_end_fasta,
    test_metaphlan_paired_end_fastq_save_sam,
]
