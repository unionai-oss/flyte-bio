"""Tests for flyte_bio.modules.trust4.

The cases and expected md5s are upstream's trust4 tests (modules commit
9a1cc9f, TRUST4 1.1.5): SARS-CoV-2 single-end and paired-end reads against
the SARS-CoV-2 transcriptome, without an IMGT reference. As upstream, the
``.out`` files and ``<prefix>_airr_align.tsv`` are only checked to exist
(upstream's unstable keys), the other outputs by md5.
"""

from flyte.io import File

from flyte_bio.modules.trust4 import Trust4Result, trust4
from tests.framework import assert_md5, env, fixture

SARSCOV2 = "genomics/sarscov2/"
TRANSCRIPTOME = SARSCOV2 + "genome/transcriptome.fasta"
AIRR_MD5 = "2cd92148076cef8ec3d475085fd4380a"
REPORT_MD5 = "51623a22f060007a69207b1ebdf3a7cd"
OUT_NAMES = ["test_cdr3.out", "test_final.out", "test_raw.out"]
TSV_NAMES = ["test_airr.tsv", "test_airr_align.tsv", "test_report.tsv"]


def name(f: File) -> str:
    return f.path.rstrip("/").rsplit("/", 1)[-1]


async def check(result: Trust4Result, fasta: dict[str, str], fq: dict[str, str]) -> None:
    await assert_md5(result.airr, AIRR_MD5, label="test_airr.tsv")
    await assert_md5(result.report, REPORT_MD5, label="test_report.tsv")
    assert [name(f) for f in result.tsv] == TSV_NAMES, f"tsv outputs: {[name(f) for f in result.tsv]}"
    assert [name(f) for f in result.out] == OUT_NAMES, f"out outputs: {[name(f) for f in result.out]}"
    for files, expected, kind in ((result.fasta, fasta, "fasta"), (result.fq, fq, "fq")):
        assert [name(f) for f in files] == list(expected), f"{kind} outputs: {[name(f) for f in files]}"
        for f in files:
            await assert_md5(f, expected[name(f)], label=name(f))


@env.task
async def test_trust4_single_end() -> None:
    reads = [await fixture(SARSCOV2 + "illumina/fastq/test_1.fastq.gz")]
    result = await trust4(await fixture(TRANSCRIPTOME), prefix="test", reads=reads, single_end=True)
    await check(
        result,
        fasta={
            "test_annot.fa": "4b07aefcac3763f7f3d242d9371060df",
            "test_assembled_reads.fa": "2c9601da894105f2fe6c068659214db3",
        },
        fq={"test_toassemble.fq": "bd7ab2f9ea96290e58e868e4aabccf1f"},
    )


@env.task
async def test_trust4_paired_end() -> None:
    reads = [
        await fixture(SARSCOV2 + "illumina/fastq/test2_1.fastq.gz"),
        await fixture(SARSCOV2 + "illumina/fastq/test2_2.fastq.gz"),
    ]
    result = await trust4(await fixture(TRANSCRIPTOME), prefix="test", reads=reads)
    await check(
        result,
        fasta={
            "test_annot.fa": "110678ccbb8f348b25396f364a777377",
            "test_assembled_reads.fa": "dce87dfb35c6efc260440865996e6c9b",
        },
        fq={
            "test_toassemble_1.fq": "bd7ab2f9ea96290e58e868e4aabccf1f",
            "test_toassemble_2.fq": "0211588305d478e10879fb35b93e13d5",
        },
    )


tests = [test_trust4_single_end, test_trust4_paired_end]
