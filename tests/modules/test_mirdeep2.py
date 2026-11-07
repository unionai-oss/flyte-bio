"""Tests for flyte_bio.modules.mirdeep2.

The cases are upstream's mirdeep2/mapper and mirdeep2/mirdeep2 tests (modules
commit 9a1cc9f, miRDeep2 2.0.1.2), with upstream's mapper ``ext.args``.
Upstream snapshots no md5s for these outputs (the collapsed read IDs and
miRDeep2's date-stamped results vary), only the file names and content checks,
which are ported as they are. Upstream's ``seqkit fq2fa`` setup step is done
here in Python; its output only feeds ``seqkit replace``, which rewrites it.
"""

import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bowtie import bowtie1_build
from flyte_bio.modules.mirdeep2 import Mirdeep2MapperResult, mirdeep2_mapper, mirdeep2_mirdeep2
from flyte_bio.modules.seqkit import seqkit_replace
from tests.framework import env, fixture

TUTORIAL = "https://github.com/rajewsky-lab/mirdeep2/raw/master/tutorial_dir/"
SMRNASEQ_READS = "https://github.com/nf-core/test-datasets/raw/smrnaseq/testdata/trimmed/small_Clone1_N1.fastp.fastq.gz"
MAPPER_ARGS = "-c -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -v"
CELEGANS_HIT = (
    "21\t1\t21\ttcaccgggtgtaaatcagctt\tchrII:11534525-11540624\t21\t3535\t3555"
    "\ttcaccgggtgtaaatcagctt\t+\t0\tmmmmmmmmmmmmmmmmmmmmm"
)


async def lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


def count(lines_: list[str], text: str) -> int:
    return sum(text in line for line in lines_)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def cel_cluster_index():
    return await bowtie1_build(await fixture(TUTORIAL + "cel_cluster.fa"), prefix="genome_cel_cluster")


async def celegans_mapping() -> Mirdeep2MapperResult:
    reads = await fixture(TUTORIAL + "reads.fa")
    return await mirdeep2_mapper(reads, await cel_cluster_index(), prefix="test_reads", args=MAPPER_ARGS)


async def fastq_to_fasta(fastq: File, name: str) -> File:
    """``seqkit fq2fa``: each record's header line and sequence, as FASTA."""
    local = await fastq.download()
    out = Path(tempfile.mkdtemp()) / name
    with gzip.open(local, "rt") as fin, out.open("w") as fout:
        for i, line in enumerate(fin):
            if i % 4 == 0:
                fout.write(">" + line[1:])
            elif i % 4 == 1:
                fout.write(line)
    return await File.from_local(str(out))


@env.task
async def test_mirdeep2_mapper_fasta_celegans() -> None:
    result = await celegans_mapping()
    assert file_name(result.collapsed) == "test_reads_collapsed.fa", result.collapsed.path
    assert file_name(result.arf) == "test_reads_reads_collapsed_vs_genome_cel_cluster_genome.arf", result.arf.path
    collapsed = await lines(result.collapsed)
    assert count(collapsed, "TCACCGGGGGTACATCAGCTAA") == 1, "TCACCGGGGGTACATCAGCTAA not collapsed to one read"
    assert count(collapsed, "seq_347479_x287") == 1, "seq_347479_x287 not found once"
    assert count(await lines(result.arf), CELEGANS_HIT) == 4, "expected 4 hits of the chrII:11534525 alignment"


@env.task
async def test_mirdeep2_mapper_fasta_smrnaseq() -> None:
    fasta = await fastq_to_fasta(await fixture(SMRNASEQ_READS), "small_Clone1_N1.fa")
    reads = await seqkit_replace(fasta, prefix="small_Clone1_N1", out_ext="fasta", args=r"-p \s.+")
    result = await mirdeep2_mapper(reads, await cel_cluster_index(), prefix="small_Clone1_N1", args=MAPPER_ARGS)
    assert file_name(result.collapsed) == "small_Clone1_N1_collapsed.fa", result.collapsed.path
    expected_arf = "small_Clone1_N1_reads_collapsed_vs_genome_cel_cluster_genome.arf"
    assert file_name(result.arf) == expected_arf, result.arf.path
    collapsed = await lines(result.collapsed)
    assert count(collapsed, "TACCTGAGGTAGCAGGTTGTATAGTTGGGG") == 1, "read not collapsed to one"
    assert count(collapsed, "seq_996152_x1") == 1, "seq_996152_x1 not found once"


@env.task
async def test_mirdeep2_mirdeep2_fa() -> None:
    mapping = await celegans_mapping()
    result = await mirdeep2_mirdeep2(
        mapping.collapsed,
        mapping.arf,
        await fixture(TUTORIAL + "cel_cluster.fa"),
        mature=await fixture(TUTORIAL + "mature_ref_this_species.fa"),
        hairpin=await fixture(TUTORIAL + "precursors_ref_this_species.fa"),
        prefix="test_reads",
    )
    names = sorted(file_name(f) for f in (result.bed, result.csv, result.html))
    assert names == ["result_test_reads.bed", "result_test_reads.csv", "result_test_reads.html"], names
    assert "</html>" in (await lines(result.html))[-1], "html report is incomplete"
    assert "miRDeep2 score" in (await lines(result.csv))[0], "csv lacks the miRDeep2 score header"


tests = [
    test_mirdeep2_mapper_fasta_celegans,
    test_mirdeep2_mapper_fasta_smrnaseq,
    test_mirdeep2_mirdeep2_fa,
]
